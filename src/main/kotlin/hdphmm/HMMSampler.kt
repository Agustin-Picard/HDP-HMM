package hdphmm

import distribution.Distribution
import distribution.GaussianDistribution1D
import tomasvolker.kyplot.dsl.showLine
import tomasvolker.kyscript.KyScriptConfig
import tomasvolker.numeriko.core.dsl.D
import tomasvolker.numeriko.core.functions.cumulativeSum
import tomasvolker.numeriko.core.interfaces.array1d.double.DoubleArray1D
import tomasvolker.numeriko.core.interfaces.factory.doubleArray1D
import tomasvolker.numeriko.core.interfaces.factory.toDoubleArray1D
import tomasvolker.numeriko.core.linearalgebra.linearSpace
import tomasvolker.numeriko.core.operations.concatenate
import tomasvolker.openrndr.math.primitives.d
import utils.autocorrelation
import utils.autocorrelationTime
import kotlin.math.absoluteValue
import kotlin.math.min
import kotlin.random.Random
import kotlin.system.measureTimeMillis

interface Sampler<P> {
    val stateCount: Int
    val learnedParameters: P
    fun initTrain(nRepetitions: Int, nIterations: Int)
    fun fit(nIterations: Int)
    fun trainStep()
}

operator fun <P,E> MarkovChain<P,out E>.get(i: Int) = getState(i)

fun ignore(): Nothing = error("illegal path")

open class HMMSampler<P,E>(val observations: DoubleArray1D,
                           private val observationDistribution: Distribution<P, E>,
                           private val priorDistribution: Distribution<*, P>,
                           val alpha: Double
): Sampler<P> {

    val duration: Int = observations.size
    val markovChain = MarkovChain<P,E>(duration).apply { initializeChain() }
    override val learnedParameters: P get() = markovChain.parameters
    override val stateCount get() = markovChain.nStates

    private fun MarkovChain<P,E>.chainDistribution(): DoubleArray1D =
            doubleArray1D(nStates + 1) {
                if (it < nStates) this[it].multiplicity / (nStateCount + alpha) else alpha / (nStateCount + alpha)
            }

    private fun MarkovChain<P,E>.initializeChain() {
        for (i in 0 until duration) {
            when (i) {
                0 -> add(proposeStateFromPrior(), i)
                else -> add(proposeState(chainDistribution()) ?: proposeStateFromPrior(), i)
            }
        }
    }

    fun MarkovChain<P,E>.proposeStateFromPrior() =
            MarkovChainState(observationDistribution.fromPrior(priorDistribution))

    private fun initTrainStep(obsPosition: Int, repetitions: Int): Unit {
        val currState = markovChain[obsPosition]

        repeat(repetitions) {
            val newState = markovChain.proposeState(markovChain.chainDistribution()) ?:
            markovChain.proposeStateFromPrior()

            if (Random.nextDouble() < min(1.0, acceptanceProbability(
                            currDist = currState.emissionDistribution,
                            newDist = newState.emissionDistribution,
                            observation = observations[obsPosition]))) {
                markovChain.updateState(newState, obsPosition)
            } else
                markovChain.updateState(currState, obsPosition)
        }
    }

    override fun initTrain(nRepetitions: Int, nIterations: Int) {
        repeat(nIterations) {
            for(i in 0 until duration) { initTrainStep(i, nRepetitions) }
        }
    }

    override fun fit(nIterations: Int) =
            repeat(nIterations) { trainStep() }

    override fun trainStep() =
        (0 until duration).forEach { initTrainStep(it, 4) }

    private fun acceptanceProbability(currDist: Distribution<P, E>, newDist: Distribution<P, E>, observation: Double) =
            newDist.probability(observation) / (currDist.probability(observation) + 1e-12)

    fun clearChain() {
        markovChain.clear()
        markovChain.initializeChain()
    }

    override fun toString(): String = markovChain.toString()
}


class NoGapsSampler<P,E>(
    observations: DoubleArray1D,
    observationDistribution: Distribution<P, E>,
    priorDistribution: Distribution<*, P>,
    alpha: Double
): HMMSampler<P,E>(
        observations,
        observationDistribution,
        priorDistribution,
        alpha
) {

    private fun MarkovChain<P,E>.chainProb(
            distinctStates: List<MarkovChainState<P,E>>,
            observation: Double
    ): DoubleArray1D =
            doubleArray1D(nStates + 1) { i ->
                if (i < nStates) distinctStates[i].multiplicity * distinctStates[i].probability(observation)
                else
                    alpha / (nStates + 1) * distinctStates.last().probability(observation)
            }.let{ it / it.sum() }

    override fun trainStep() =
        (0 until duration).forEach { i ->
            if (Random.nextDouble() > markovChain.nStates / (markovChain.nStates + 1.0)) {
                val currState = markovChain[i]

                val distinctStates = markovChain.stateList +
                        if (currState.isSingleton) currState
                        else markovChain.proposeStateFromPrior()

                markovChain.proposeState(
                        chainDistribution = markovChain.chainProb(distinctStates, observations[i]),
                        state = distinctStates.last()
                ).apply { markovChain.updateState(this!!, i) }
            }
        }
}


class AuxGibbsSampler<P,E>(
    observations: DoubleArray1D,
    observationDistribution: Distribution<P, E>,
    priorDistribution: Distribution<*, P>,
    alpha: Double,
    private val nAuxiliarStates: Int
): HMMSampler<P,E>(
        observations,
        observationDistribution,
        priorDistribution,
        alpha
) {
    private fun MarkovChain<P,E>.chainProbability(
            states: List<MarkovChainState<P,E>>,
            observation: Double): DoubleArray1D =
        doubleArray1D(states.size) { i ->
            if (i < nStates - nAuxiliarStates)
                states[i].multiplicity.toDouble() * states[i].probability(observation)
            else
                (alpha / nAuxiliarStates) / (nStates + alpha - 1.0) * states[i].probability(observation)
        }.let{ it / it.sum() }

    private fun MarkovChain<P,E>.proposeStateFromAuxChain(
            states: List<MarkovChainState<P, E>>,
            observation: Double): MarkovChainState<P, E> {
        val randU = Random.nextDouble()
        val index = List(states.size) { i ->
            states[i] to chainProbability(states, observation).cumulativeSum()[i]
        }.asSequence()
                .sortedBy { it.second }
                .indexOfFirst { it.second > randU }

        return if (index == -1)
            states.last()
        else
            states[index]
    }

    private fun MarkovChain<P,E>.proposeStateFromPrior(nStates: Int) =
            List(nStates) { proposeStateFromPrior() }

    override fun trainStep() {
        for (i in 0 until duration) {
            markovChain.updateState(
                    state = markovChain.proposeStateFromAuxChain(
                            states = markovChain.stateList + proposeAuxStates(markovChain[i]),
                            observation = observations[i]),
                    observation = i
            )
        }
    }

    private fun proposeAuxStates(currState: MarkovChainState<P,E>) =
            if (currState.isSingleton)
                listOf(currState) + markovChain.proposeStateFromPrior(nAuxiliarStates - 1)
            else
                markovChain.proposeStateFromPrior(nAuxiliarStates)
}


class PartialGibbsHMSampler<P,E>(
    observations: DoubleArray1D,
    observationDistribution: Distribution<P, E>,
    priorDistribution: Distribution<*, P>,
    alpha: Double
): HMMSampler<P,E>(
        observations,
        observationDistribution,
        priorDistribution,
        alpha
) {
    private fun MarkovChain<P,E>.chainProbWeighted(obsPosition: Int) =
            doubleArray1D(nStates) { i ->
                if (stateList[i].contains(obsPosition))
                    (stateList[i].probability(observations[obsPosition] * stateList[i].multiplicity) / (duration - 1.0))
                else
                    (stateList[i].probability(observations[obsPosition]) * (stateList[i].multiplicity - 1)) / (duration - 1.0)
            }.let{ it / it.sum() }

    private fun MarkovChain<P,E>.chainProb() =
            doubleArray1D(nStates) { i ->
                stateList[i].multiplicity / nStateCount
            }

    private fun acceptanceProbability(currState: MarkovChainState<P,E>, newState: MarkovChainState<P,E>, observation: Double) =
            newState.probability(observation) / currState.probability(observation)

    override fun trainStep() {
        for (i in 0 until duration) {
            val currState = markovChain[i]
            if (markovChain[i].isSingleton) {
                val newState = markovChain.proposeState(markovChain.chainProb()) ?: markovChain.stateList.last()

                if (min(1.0, acceptanceProbability(currState, newState, observations[i]) * alpha / (stateCount - 1.0)) > Random.nextDouble())
                    markovChain.updateState(newState, i)
            } else {
                val newState = markovChain.proposeStateFromPrior()

                if (acceptanceProbability(currState, newState, observations[i]) * (stateCount - 1.0) / alpha > Random.nextDouble())
                    markovChain.updateState(newState, i)
            }

            if (markovChain[i].isSingleton)
                with(markovChain) { updateState(proposeState(chainProbWeighted(i)) ?: ignore(), i) }
        }
    }
}


class HMMSamplerTester<P>(private val sampler: Sampler<P>, nInitIterations: Int) {

    init {
        sampler.initTrain(nRepetitions = 5, nIterations = nInitIterations)
        println("Initialized chain: $sampler")
    }

    fun testParameter(nIterations: Int): List<P> {
        val parameterList = mutableListOf<P>()

        for (i in 0 until nIterations) {
            sampler.trainStep()
            parameterList.add(sampler.learnedParameters)
        }

        return parameterList
    }

    fun testStateCount(nIterations: Int): List<Int> {
        val stateCountList = mutableListOf<Int>()

        for (i in 0 until nIterations) {
            sampler.trainStep()
            stateCountList.add(sampler.stateCount)
        }

        return stateCountList
    }

    fun testMeasureTime(nIterations: Int): Double {
        val iterationTimeList = mutableListOf<Long>()

        for (i in 0 until nIterations) {
            iterationTimeList.add(measureTimeMillis { sampler.trainStep() })
        }

        return iterationTimeList.average().also { println(it) }
    }

    fun plotAutocorrelation(samples: DoubleArray1D, start: Int = 0, stop: Int = samples.size) {
        val newArray = samples[0 until start].average().let {
            listOf(it).toDoubleArray1D()
                .concatenate(samples[start until stop])
        }
        val autocorrList = List(stop - start) { i -> newArray.autocorrelation(i) }

        showLine(title = "Autocorrelation") {
            x = linearSpace(start.d, stop.d, stop - start)
            y = autocorrList
        }
    }

    override fun toString(): String = sampler.toString()
}


fun main() {
    KyScriptConfig.defaultPythonPath = "python"

    val observations = D[-1.48, -1.40, -1.16, -1.08, -1.02, 0.14, 0.51, 0.53, 0.78]
    val observationDistribution = GaussianDistribution1D(mean = 0.0, std = 0.1)
    val priorDistribution = GaussianDistribution1D(mean = 0.0, std = 1.0)

    /*val noGapsSampler = NoGapsSampler(
            observations = observations,
            observationDistribution = observationDistribution,
            priorDistribution = priorDistribution,
            alpha = 1.0
    )*/

    val auxGibbsSampler = AuxGibbsSampler(
            observations = observations,
            observationDistribution = observationDistribution,
            priorDistribution = priorDistribution,
            alpha = 1.0,
            nAuxiliarStates = 30
    )

    /*val partialGibbsHMSampler = PartialGibbsHMSampler(
            observations = observations,
            observationDistribution = observationDistribution,
            priorDistribution = priorDistribution,
            alpha = 1.0
    )*/

    /*val metropolisHastingsSampler = HMMSampler(observations = observations,
            observationDistribution = observationDistribution,
            priorDistribution = priorDistribution,
            alpha = 1.0
    )*/

    /*val tester: HMMSamplerTester<Double> = HMMSamplerTester(auxGibbsSampler, 100).apply {
        val parameter = testParameter(40000).toDoubleArray1D()
        println("Autocorrelation time = ${parameter.autocorrelationTime(30000)}")
        println(this)
        plotAutocorrelation(parameter, 30000)
    }*/

    val tester: HMMSamplerTester<Double> = HMMSamplerTester(auxGibbsSampler, 100).apply {
        val autocorrList = mutableListOf<Double>()

        repeat(1000) {
//            val parameter = testMeasureTime(20000)
            val parameter = testStateCount(20000).toDoubleArray1D()
            autocorrList.add(parameter.autocorrelationTime(15000).absoluteValue)
            println(it)
        }

        println(autocorrList.average())
    }

}