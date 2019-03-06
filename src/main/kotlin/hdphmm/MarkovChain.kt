package hdphmm

import distribution.Distribution
import tomasvolker.numeriko.core.functions.cumulativeSum
import tomasvolker.numeriko.core.interfaces.array1d.double.DoubleArray1D
import kotlin.random.Random

data class MarkovChainState<P,out E>(val emissionDistribution: Distribution<P, E>) {

    lateinit var observations: MutableList<Int>

    constructor(emissionDistribution: Distribution<P, E>,
                observation: Int): this(emissionDistribution) {
        if (!::observations.isInitialized)
            observations = mutableListOf(observation)
        else {
            observations.add(observation)
        }
    }

    constructor(emissionDistribution: Distribution<P, E>,
                observation: List<Int>): this(emissionDistribution) {
        if (!::observations.isInitialized)
            observations = observation.toMutableList()
        else {
            this.observations.addAll(observations)
        }
    }

    val multiplicity: Int
        get() = if (!::observations.isInitialized) 0 else observations.size
    val isSingleton get() = multiplicity == 1
    val parameters get() = emissionDistribution.parameters

    fun newValue(): E = emissionDistribution.nextSample()

    fun addObservation(observation: Int) {
        if (!::observations.isInitialized) observations = mutableListOf(observation)
        else {
            if (observation !in observations) observations.add(observation)
        }
    }

    fun addObservation(observation: List<Int>) {
        if (!::observations.isInitialized) observations = observation.toMutableList()
        else {
            this.observations.addAll(observations.filter { it !in this.observations})
        }
    }

    fun removeObservation(observation: Int) {
        if (!::observations.isInitialized) observations = mutableListOf(observation)
        else {
            if (observation in observations) observations.remove(element = observation)
        }
    }

    fun removeObservation(observation: List<Int>) {
        if (!::observations.isInitialized) observations = observation.toMutableList()
        else {
            this.observations.filter { it in observation }
        }
    }

    fun probability(observation: Double) =
            emissionDistribution.probability(observation)

    fun contains(observation: Int) =
            observations.contains(observation)

    override fun toString(): String = "Distribution with params $parameters and observations $observations"
}

class MarkovChain<P,E>(val chainSize: Int) {

    val stateList: MutableList<MarkovChainState<P,E>> = mutableListOf()
    val nStates: Int get() = stateList.size
    val nStateCount: Int get() = countStates()
    val parameters: P get() = stateList[0].parameters

    fun clear() = stateList.clear()

    fun add(state: MarkovChainState<P, E>, observation: Int) {
        val newState = getState(state)
        if (newState == null) {
            state.addObservation(observation)
            stateList.add(state)
        } else {
            newState.addObservation(observation)
        }
    }

    fun remove(observation: Int): Unit {
        val state = getState(observation)
        if (state.observations.size <= 1)
            stateList.remove(state)
        else
            state.removeObservation(observation)
    }

    fun updateState(state: MarkovChainState<P,E>, observation: Int) {
        remove(observation)
        state.addObservation(observation)
        add(state, observation)
    }

    fun contains(state: MarkovChainState<P,E>): Boolean {
        return stateList.filter { it.emissionDistribution == state.emissionDistribution }.isNotEmpty()
    }

    fun getState(state: MarkovChainState<P,E>): MarkovChainState<P,E>? =
            stateList.find { it.emissionDistribution == state.emissionDistribution }

    fun getState(observation: Int): MarkovChainState<P,E> =
            stateList.find { it.observations.contains(observation) } ?: error("observation not contained")

    // Algorithm taken from: https://stats.stackexchange.com/questions/67911/how-to-sample-from-a-discrete-distribution
    fun proposeState(chainDistribution: DoubleArray1D,
                     random: Random = Random.Default,
                     state: MarkovChainState<P, E>? = null): MarkovChainState<P,E>? {
        val randU = random.nextDouble()
        val chainCumSum = chainDistribution.cumulativeSum()
        val index =  List(chainDistribution.size) { i->
            if (i < nStates)
                stateList[i] to chainCumSum[i]
            else
                state to chainCumSum[i]
        }.asSequence()
                .sortedBy { it.second }
                .indexOfFirst { it.second > randU }

        return if (index > nStates - 1 || index == -1)
            state
        else
            stateList[index]
    }

    fun countStates(): Int {
        var result: Int = 0

        stateList.forEach { result += it.observations.count() }

        return result
    }

    override fun toString(): String {
        var result: String = ""
        stateList.forEach { result += it.toString() + "\n" }

        return result
    }

}