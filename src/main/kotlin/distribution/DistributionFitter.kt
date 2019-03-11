package distribution

import org.apache.commons.math3.analysis.UnivariateFunction
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction
import org.apache.commons.math3.analysis.solvers.*
import tomasvolker.numeriko.core.functions.average
import utils.mean
import tomasvolker.numeriko.core.interfaces.array1d.double.DoubleArray1D
import tomasvolker.numeriko.core.interfaces.array1d.double.elementWise
import tomasvolker.numeriko.core.interfaces.factory.toDoubleArray1D
import tomasvolker.numeriko.core.primitives.squared
import tomasvolker.numeriko.core.primitives.sumDouble
import utils.digammaFunction
import utils.trigammaFunction
import kotlin.math.absoluteValue
import kotlin.math.ln
import kotlin.math.pow
import kotlin.math.sqrt

fun DoubleArray1D.std(): Double =
        mean().let { mean ->
            sqrt(sumDouble(0 until size) { i ->
                (this[i] - mean).squared() } / (size - 1)) }

interface DistributionFitter<E,P> {
    fun fitParameters(samples: Iterable<E>): P
}

data class GaussianDistribution1DParameters(val mean: Double, val std: Double)

object GaussianDistributionFitter1D: DistributionFitter<Double, GaussianDistribution1DParameters> {
    override fun fitParameters(samples: Iterable<Double>): GaussianDistribution1DParameters =
            samples.toDoubleArray1D().let {
                GaussianDistribution1DParameters(
                    mean = estimateMean(it),
                    std = estimateStd(it)
                )
            }

    private fun estimateMean(samples: DoubleArray1D) =
            samples.average()

    private fun estimateStd(samples: DoubleArray1D) =
            samples.std()
}

data class PoissonDistributionParameters(val rate: Double)

object PoissonDistributionFitter: DistributionFitter<Double, PoissonDistributionParameters> {
    private fun estimateRate(samples: DoubleArray1D) =
            samples.sum() / samples.size

    override fun fitParameters(samples: Iterable<Double>): PoissonDistributionParameters =
        PoissonDistributionParameters(
            estimateRate(
                samples.toDoubleArray1D()
            )
        )
}

data class LaplaceDistributionParameters(val mean: Double, val scale: Double)

object LaplaceDistributionFitter: DistributionFitter<Double, LaplaceDistributionParameters> {
    override fun fitParameters(samples: Iterable<Double>): LaplaceDistributionParameters =
            samples.toDoubleArray1D().let {
                LaplaceDistributionParameters(
                    mean = estimateMean(it),
                    scale = estimateScale(it)
                )
            }

    private fun estimateMean(samples: DoubleArray1D) =
            samples.average()

    private fun estimateScale(samples: DoubleArray1D) =
            2 * samples.std().squared()
}

data class BetaDistributionParameters(val alpha: Double, val beta: Double)

object BetaDistributionFitter: DistributionFitter<Double, BetaDistributionParameters> {
    override fun fitParameters(samples: Iterable<Double>): BetaDistributionParameters {
        val mean = samples.average()
        val std = samples.toDoubleArray1D().std()

        return BetaDistributionParameters(
            alpha = estimateAlpha(
                mean,
                std
            ), beta = estimateBeta(mean, std)
        )
    }

    private fun estimateAlpha(mean: Double, std: Double) =
            if (std.squared() < mean * (1 - mean))
                mean * (mean * (1.0 - mean) / std.squared() - 1.0)
            else
                0.0

    private fun estimateBeta(mean: Double, std: Double) =
            if (std.squared() < mean * (1 - mean))
                (1.0 - mean) * ((mean * (1.0 - mean)) / std.squared() - 1.0)
            else
                0.0
}

data class GeometricDistributionParameters(val successProbability: Double)

object GeometricDistributionFitter: DistributionFitter<Double, GeometricDistributionParameters> {
    override fun fitParameters(samples: Iterable<Double>): GeometricDistributionParameters =
        GeometricDistributionParameters(
            estimateSuccessProbability(
                samples.toDoubleArray1D()
            )
        )

    private fun estimateSuccessProbability(samples: DoubleArray1D) =
            1.0 / samples.mean()
}

data class ExponentialDistributionParameters(val decayRate: Double)

object ExponentialDistributionFitter: DistributionFitter<Double, ExponentialDistributionParameters> {
    override fun fitParameters(samples: Iterable<Double>) =
        ExponentialDistributionParameters(
            estimateDecayRate(
                samples.toDoubleArray1D()
            )
        )

    private fun estimateDecayRate(samples: DoubleArray1D) =
            (samples.size - 2.0) / samples.sum()

}

data class GammaDistributionParameters(val scale: Double, val shape: Double)

object GammaDistributionFitter: DistributionFitter<Double, GammaDistributionParameters> {

    override fun fitParameters(samples: Iterable<Double>): GammaDistributionParameters {
        val samplesMean = samples.average()
        val samplesStd = samples.toDoubleArray1D().std()
        val shape = (samplesMean / samplesStd).squared()
        val scale = samplesStd.squared() / samplesMean

        return GammaDistributionParameters(scale = scale, shape = shape)
    }

    fun fitParametersMLE(samples: Iterable<Double>): GammaDistributionParameters {
        val meanLog = samples.map { ln(it) }.average()
        val logMean = ln(samples.average())
        var shape = 1.0
        var currShape = 0.0

        while ((shape - currShape).absoluteValue > 1e-2) {
            currShape = shape
            shape = 1.0 / ((1 / shape) + (meanLog - logMean + ln(shape) - digammaFunction(shape)) /
                    (shape.squared() * ((1 / shape) - trigammaFunction(shape))))
        }

        return GammaDistributionParameters(scale = samples.average() / shape, shape = shape)
    }

    fun fitParametersMoments(samples: Iterable<Double>): GammaDistributionParameters {
        val samplesMean = samples.average()
        val samplesStd = samples.toDoubleArray1D().std()
        val shape = (samplesMean / samplesStd).squared()
        val scale = samplesStd.squared() / samplesMean

        return GammaDistributionParameters(scale = scale, shape = shape)
    }

}

data class ParetoDistributionParameters(val scale: Double, val alpha: Double)

object ParetoDistributionFitter: DistributionFitter<Double, ParetoDistributionParameters> {
        override fun fitParameters(samples: Iterable<Double>): ParetoDistributionParameters =
                samples.toDoubleArray1D().let {
                    ParetoDistributionParameters(
                        scale = estimateScale(it),
                        alpha = estimateAlpha(
                            it,
                            estimateScale(it)
                        )
                    )
                }

        private fun estimateScale(samples: DoubleArray1D) = samples.min() ?: 0.0

        private fun estimateAlpha(samples: DoubleArray1D, scale: Double) =
                samples.size / sumDouble(0 until samples.size) { i -> ln(samples[i] / scale) }
    }

data class LogNormalDistributionParameters(val mean: Double, val std: Double)

object LogNormalDistributionFitter: DistributionFitter<Double, LogNormalDistributionParameters> {
        override fun fitParameters(samples: Iterable<Double>) =
            samples.toDoubleArray1D().let {
                LogNormalDistributionParameters(
                    mean = estimateMean(it),
                    std = estimateVariance(it)
                )
            }

        private fun estimateMean(samples: DoubleArray1D) = samples.elementWise { ln(it) }.mean()

        private fun estimateVariance(samples: DoubleArray1D) =
                    samples.elementWise { ln(it) }.std()
    }


data class WeibullDistributionParameters(val scale: Double, val shape: Double)

object WeibullDistributionFitter: DistributionFitter<Double, WeibullDistributionParameters> {
    override fun fitParameters(samples: Iterable<Double>): WeibullDistributionParameters =
            samples.toDoubleArray1D().let {
                val scale = estimateScale(it)
                WeibullDistributionParameters(
                    scale = scale,
                    shape = estimateShape(it, scale)
                )
            }

    private fun estimateShape(samples: DoubleArray1D, scale: Double): Double =
        (samples.elementWise { it.pow(scale) }.sum() / samples.size).pow(1 / scale)

    private fun estimateScale(samples: DoubleArray1D): Double {
        val shapeFunction: UnivariateFunction = UnivariateFunction {
                k: Double -> sumDouble(0 until samples.size) { i -> samples[i].pow(k) * ln(samples[i]) } /
                sumDouble(0 until samples.size) { i -> samples[i].pow(k) } - (1.0 / k)
            - samples.elementWise { ln(it) }.average()
        }

        return BisectionSolver().solve(1000, shapeFunction, 0.0, 100.0) // Brent or Secant solver don't seem to work
    }

}