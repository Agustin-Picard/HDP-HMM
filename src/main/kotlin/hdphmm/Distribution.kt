package hdphmm

import tomasvolker.numeriko.core.interfaces.array1d.double.DoubleArray1D
import tomasvolker.numeriko.core.interfaces.arraynd.double.DoubleArrayND
import tomasvolker.numeriko.core.interfaces.factory.doubleArray1D
import tomasvolker.numeriko.core.interfaces.factory.nextGaussian
import tomasvolker.numeriko.core.primitives.squared
import tomasvolker.numeriko.core.primitives.sumDouble
import kotlin.math.*
import kotlin.random.Random

interface Distribution<P,out E> {
    val parameters: P
    fun nextSample(random: Random = Random.Default): E
    fun probability(observation: Double): Double
    fun fromPrior(distribution: Distribution<*,P>, random: Random = Random.Default): Distribution<P,E>
    fun pdf(nSamples: Int, start: Double, stop: Double): DoubleArrayND
    fun quantile(probability: Double): Double
}

interface DoubleDistribution<P>: Distribution<P,Double>{
    override fun nextSample(random: Random): Double = nextDouble(random)
    fun nextDouble(random: Random = Random.Default): Double
    override fun pdf(nSamples: Int, start: Double, stop: Double): DoubleArray1D
}

// KLD = - sum(p * ln(q / p))
fun DoubleDistribution<*>.estimateKLDivergence(other: DoubleDistribution<*>, resolution: Int = 1000): Double {
    val startPdf = min(quantile(0.01), other.quantile(0.01))
    val stopPdf = max(quantile(0.99), other.quantile(0.99))

    val thisPdf = pdf(resolution, startPdf, stopPdf)
    val otherPdf = other.pdf(resolution, startPdf, stopPdf)

    return sumDouble(0 until resolution) { i -> thisPdf[i] * ln(otherPdf[i] / thisPdf[i]) }
}

fun DoubleDistribution<*>.estimateEntropy(resolution: Int = 1000): Double =
        pdf(resolution, quantile(0.01), quantile(0.99)).let { pdf ->
            sumDouble(0 until resolution) { i ->  pdf[i] * ln(pdf[i]) }
        }

// Classical 1D normal distribution with input parameter mean and 1D output
class GaussianDistribution(val mean: Double, val std: Double): DoubleDistribution<Double> {

    override val parameters: Double get() = mean
    val variance: Double get() = std.squared()

    override fun nextDouble(random: Random): Double =
        mean + std * random.nextGaussian()

    fun nextSample(nSamples: Int, random: Random = Random.Default): DoubleArray1D =
            doubleArray1D(nSamples) { nextDouble(random) }

    override fun probability(observation: Double): Double =
            exp(- (observation - mean).squared() / (2 * variance)) / (sqrt(2 * PI) * std)

    override fun fromPrior(distribution: Distribution<*, Double>, random: Random): Distribution<Double, Double> =
            GaussianDistribution(
                    mean = distribution.nextSample(random) + mean,
                    std = std)

    override fun pdf(nSamples: Int, start: Double, stop: Double): DoubleArray1D =
        doubleArray1D(nSamples) {
            exp(-((it * (stop - start) / nSamples) - mean).squared() / (2 * variance)) / (sqrt(2 * PI) * std)
        }

    override fun quantile(probability: Double): Double =
            TODO("calculate the erfinv...")
}