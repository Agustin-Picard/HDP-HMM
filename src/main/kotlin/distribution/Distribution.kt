package distribution

import org.openrndr.math.Vector2
import tomasvolker.numeriko.core.interfaces.array1d.double.DoubleArray1D
import tomasvolker.numeriko.core.interfaces.arraynd.double.DoubleArrayND
import tomasvolker.numeriko.core.interfaces.factory.doubleArray1D
import tomasvolker.numeriko.core.interfaces.factory.nextGaussian
import tomasvolker.numeriko.core.primitives.squared
import tomasvolker.numeriko.core.primitives.sumDouble
import tomasvolker.openrndr.math.plot.plotLine
import tomasvolker.openrndr.math.plot.quickPlot2D
import utils.erfinvFunction
import utils.gammaFunction
import kotlin.math.*
import kotlin.random.Random

interface Distribution<P,out E> {
    val parameters: P
    fun nextSample(random: Random = Random.Default): E
    fun probability(observation: Double): Double
    fun fromPrior(distribution: Distribution<*, P>, random: Random = Random.Default): Distribution<P, E>
    fun pdf(nSamples: Int, start: Double, stop: Double): DoubleArrayND
    fun quantile(probability: Double): Double
}

interface DoubleDistribution<P>: Distribution<P, Double> {
    override fun nextSample(random: Random): Double = nextDouble(random)
    fun nextSample(nSamples: Int, random: Random = Random.Default) =
            doubleArray1D(nSamples) { nextDouble(random) }
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
class GaussianDistribution1D(val mean: Double, val std: Double): DoubleDistribution<Double> {

    override val parameters: Double get() = mean
    val variance: Double get() = std.squared()

    override fun nextDouble(random: Random): Double =
        mean + std * random.nextGaussian()

    override fun probability(observation: Double): Double =
            exp(- (observation - mean).squared() / (2 * variance)) / (sqrt(2 * PI) * std)

    override fun fromPrior(distribution: Distribution<*, Double>, random: Random): Distribution<Double, Double> =
        GaussianDistribution1D(
            mean = distribution.nextSample(random) + mean,
            std = std
        )

    override fun pdf(nSamples: Int, start: Double, stop: Double): DoubleArray1D =
        doubleArray1D(nSamples) {t ->
            val x = (stop - start)  / nSamples * t + start
            probability(x)
        }

    override fun quantile(probability: Double): Double =
            mean + std * sqrt(2.0) * erfinvFunction(2.0 * probability - 1.0)

    companion object {
        fun fromParameters(param: GaussianDistribution1DParameters) =
                GaussianDistribution1D(param.mean, param.std)
    }
}

class ExponentialDistribution(val rate: Double): DoubleDistribution<Double> {

    override val parameters get() = rate
    val mean: Double get() = 1 / rate
    val variance: Double get() = rate.pow(-2)

    override fun nextDouble(random: Random): Double =
            - ln(random.nextDouble()) / rate

    override fun fromPrior(distribution: Distribution<*, Double>, random: Random): Distribution<Double, Double> =
            ExponentialDistribution(
                rate = distribution.nextSample(random)
            )

    override fun pdf(nSamples: Int, start: Double, stop: Double): DoubleArray1D =
            doubleArray1D(nSamples) { t ->
                val x = (stop - start) * t / nSamples + start
                probability(x)
            }

    override fun probability(observation: Double): Double =
            rate * exp(-rate * observation)

    override fun quantile(probability: Double): Double =
            -ln(1.0 - probability) / rate

    companion object {
        fun fromParameters(param: ExponentialDistributionParameters) =
                ExponentialDistribution(param.decayRate)
    }
}

class LogNormalDistribution(val mean: Double, val std: Double): DoubleDistribution<LogNormalDistributionParameters> {

    override val parameters: LogNormalDistributionParameters
        get() = LogNormalDistributionParameters(mean, std)
    val variance: Double get() = std.squared()

    override fun nextDouble(random: Random): Double =
            ln(random.nextGaussian())

    override fun fromPrior(
        distribution: Distribution<*, LogNormalDistributionParameters>,
        random: Random
    ): Distribution<LogNormalDistributionParameters, Double> =
        distribution.nextSample().let {
            LogNormalDistribution(it.mean, it.std)
        }

    override fun pdf(nSamples: Int, start: Double, stop: Double): DoubleArray1D =
            doubleArray1D(nSamples) { t ->
                val x = (stop - start) * t / nSamples + start
                probability(x)
            }

    override fun probability(observation: Double): Double =
            (1 / (observation * std * sqrt(2 * PI))) *
                    exp(-(ln(observation) - mean).squared() / (2.0 * std).squared())

    override fun quantile(probability: Double): Double =
            exp(mean + std * sqrt(2.0) * erfinvFunction(2.0 * probability - 1.0))

    companion object {
        fun fromParameters(param: LogNormalDistributionParameters) =
                LogNormalDistribution(param.mean, param.std)
    }
}

class GammaDistribution(val shape: Double, val scale: Double): DoubleDistribution<GammaDistributionParameters> {

    override val parameters: GammaDistributionParameters
        get() = GammaDistributionParameters(scale, shape)
    val mean get() = scale * shape
    val variance get() = shape * scale.squared()

    override fun nextDouble(random: Random): Double {
        val rest = shape % 1.0
        val delta = shape - rest
        val randomArray = doubleArray1D(rest.toInt()) { random.nextDouble() }
        var eta = 10.0
        var xi = 0.0

        while (eta > xi.pow(delta - 1.0) * exp(-xi)) {
            val u = random.nextDouble()
            val v = random.nextDouble()
            val w = random.nextDouble()

            if (u <= E / (E + delta)) {
                xi = v.pow(1 / delta)
                eta = w * xi.pow(delta - 1.0)
            } else {
                xi = 1 - ln(v)
                eta = w * exp(-xi)
            }
        }

        return scale * (xi - randomArray.sumByDouble { ln(it) })
    }

    override fun fromPrior(
        distribution: Distribution<*, GammaDistributionParameters>,
        random: Random
    ): Distribution<GammaDistributionParameters, Double> =
            distribution.nextSample(random).let {
                GammaDistribution(
                    scale = it.scale,
                    shape = it.shape
                )
            }

    override fun pdf(nSamples: Int, start: Double, stop: Double): DoubleArray1D {
        val k = (1.0 / (gammaFunction(shape) * scale.pow(shape)))

        return doubleArray1D(nSamples) { t ->
            val x = (stop - start) * t / nSamples + start
            k * x.pow(shape - 1) * exp(-x / scale)
        }
    }

    override fun probability(observation: Double): Double =
        (1.0 / (gammaFunction(shape) * scale.pow(shape))) * observation.pow(shape - 1) * exp(-observation / scale)

    override fun quantile(probability: Double): Double {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    companion object {
        fun fromParameters(param: GammaDistributionParameters) =
                GammaDistribution(param.shape, param.scale)
    }
}