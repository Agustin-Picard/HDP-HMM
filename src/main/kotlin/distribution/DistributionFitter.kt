package distribution

import hdphmm.mean
import tomasvolker.numeriko.core.interfaces.array1d.double.DoubleArray1D
import tomasvolker.numeriko.core.interfaces.factory.toDoubleArray1D
import tomasvolker.numeriko.core.primitives.squared
import tomasvolker.numeriko.core.primitives.sumDouble
import kotlin.math.ln
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
            1.0 / samples.mean()

}

data class GammaDistributionParameters(val theta: Double,
                                       val k: Double,
                                       val alpha: Double,
                                       val beta: Double)

object GammaDistributionFitter: DistributionFitter<Double, GammaDistributionParameters> {
        override fun fitParameters(samples: Iterable<Double>): GammaDistributionParameters {
            val s = estimateS(samples.toDoubleArray1D())
            val k = estimateK(s)
            val theta = estimateTheta(samples.toDoubleArray1D(), k)
            val beta = 1.0 / theta

            return GammaDistributionParameters(
                theta = theta,
                k = k,
                alpha = k,
                beta = beta
            )
        }

        private fun estimateS(samples: DoubleArray1D) =
                ln(samples.mean()) - sumDouble(0 until samples.size) { i -> ln(samples[i]) } / samples.size

        private fun estimateK(s: Double) =
                (3.0 - s + sqrt((s - 3.0).squared() + 24 * s)) / (12 * s)

        private fun estimateTheta(samples: DoubleArray1D, k: Double) =
                samples.mean() / k

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
                    std = estimateStd(it)
                )
            }

        private fun estimateMean(samples: DoubleArray1D) = samples.map { ln(it) }.sum() / samples.size

        private fun estimateStd(samples: DoubleArray1D) =
                estimateMean(samples).let { mean ->
                    samples.map { (ln(it) - mean).squared() }.sum() / samples.size
                }
    }
