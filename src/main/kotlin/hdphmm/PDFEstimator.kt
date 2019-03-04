package hdphmm

import tomasvolker.kyplot.dsl.*
import tomasvolker.kyscript.KyScriptConfig
import tomasvolker.numeriko.core.dsl.D
import tomasvolker.numeriko.core.functions.*
import tomasvolker.numeriko.core.index.All
import tomasvolker.numeriko.core.interfaces.array1d.double.DoubleArray1D
import tomasvolker.numeriko.core.interfaces.array1d.double.elementWise
import tomasvolker.numeriko.core.interfaces.array2d.double.DoubleArray2D
import tomasvolker.numeriko.core.interfaces.array2d.double.elementWise
import tomasvolker.numeriko.core.interfaces.factory.*
import tomasvolker.numeriko.core.linearalgebra.linearSpace
import tomasvolker.numeriko.core.operations.stack
import tomasvolker.numeriko.core.operations.unstack
import tomasvolker.numeriko.core.primitives.squared
import tomasvolker.numeriko.core.primitives.sumDouble
import utils.singularValueDecomposition
import kotlin.math.PI
import kotlin.math.exp
import kotlin.math.pow
import kotlin.math.sqrt
import kotlin.random.Random

typealias Kernel<T> = (T) -> T
typealias DoubleKernel = (Double) -> Double

private fun DoubleArray2D.pow(x: Double) =
        singularValueDecomposition().let {
            it.u matMul it.s.elementWise { it.pow(x) } matMul it.v.transpose()
        }

private fun DoubleArray1D.normSquared() = elementWise { it.squared() }.sum()

fun DoubleArray1D.estimateKurtosis() =
        mean().let { mean ->
            sumDouble(0 until size) { i -> (this[i] - mean).pow(4) } /
                    (sumDouble(0 until size) { i -> (this[i] - mean).squared() }.squared() / size) - 3.0
        }

fun gaussianKernel1D(x: Double) = exp(-(x * x) / 2.0) / sqrt(2.0 * PI )

fun sphericalKernel1D(x: Double) = if (x <= 1.0) 1.0 else 0.0

fun gaussianKernel2D(x: DoubleArray1D) =
        exp(x.elementWise { -it.squared() / 2.0 } ) / sqrt((2.0 * PI).pow(x.size))

fun sumVector(indices: IntProgression, selector: (Int) -> DoubleArray1D) =
        indices.asSequence()
                .map { selector(it) }
                .reduce { acc, next -> acc + next }


class DoubleKernelDensityEstimator {
    companion object {
        inline fun estimatePdf(samples: Iterable<Number>,
                        start: Double,
                        stop: Double,
                        pdfCount: Int,
                        crossinline kernel: DoubleKernel,
                        bandwidth: Double) =
                SampledDistribution1D(
                        distribution = samples.toDoubleArray1D().let {
                            doubleArray1D(pdfCount) { t ->
                                sumDouble(0 until it.size) { i ->
                                    kernel(((start + t * (stop - start) / pdfCount) - it[i]) / bandwidth)
                                }
                        } / (bandwidth * it.size) },
                        start = start,
                        stop = stop
                )

        fun estimatePdf(samples: Iterable<Number>,
                        start: Double,
                        stop: Double,
                        pdfCount: Int,
                        kernel: DoubleKernel): DoubleArray1D {
            val sampleArray = samples.toDoubleArray1D()
            val bandwidth = silvermanBandwidth(sampleArray)

            return doubleArray1D(pdfCount) { t ->
                sumDouble(0 until sampleArray.size) { i ->
                    kernel(((start + t * (stop - start) / pdfCount) - sampleArray[i]) / bandwidth)
                }
            } / (bandwidth * sampleArray.size)
        }

        private fun silvermanBandwidth(samples: DoubleArray1D) =
                (4.0 * samples.std() / (3.0 * samples.size)).pow(0.2)
    }
}


class KernelDensityEstimator {
    companion object {
        fun estimatePdf(samples: DoubleArray2D,
                        start: DoubleArray1D,
                        stop: DoubleArray1D,
                        pdfCount: Int,
                        kernel: Kernel<DoubleArray1D>): SampledDistribution2D {
            val invBandwidth = silvermanBandwidth(samples).inverse()
            val detInvBandwidth = invBandwidth.determinant()

            return SampledDistribution2D(
                    distribution = List(pdfCount) { t ->
                        sumVector(0 until samples.shape1) { i ->
                            kernel(invBandwidth matMul
                                    (start + (stop - start).elementWise { t * it / pdfCount } - samples[All, i]))
                        }
                    }.map { it * detInvBandwidth }
                            .stack(),
                    start = start,
                    stop = stop
            )
        }

        private fun silvermanBandwidth(samples: DoubleArray2D) =
                samples.unstack().map { it.std() }.let { stdVector ->
                    doubleArray2D(samples.shape0, samples.shape0) {  i0, i1 ->
                        if (i0 == i1)
                            (4.0 / (samples.shape0 + 2.0)).pow(1.0 / (samples.shape0 + 4.0)) *
                                    (samples.shape1.toDouble()).pow(-1.0 / (samples.shape0 + 4.0)) *
                                    stdVector[i0]
                        else
                            0.0
                    }
                }

        private fun scottBandwidth(samples: DoubleArray2D) =
                samples.unstack().map { it.std() }.let { stdVector ->
                    doubleArray2D(samples.shape0, samples.shape0) { i0, i1 ->
                        if (i0 == i1)
                            samples.shape1.toDouble().pow(-1.0 / (samples.shape0 + 4)) * stdVector[i0]
                        else
                            0.0
                    }
                }
    }
}


fun main() {
    KyScriptConfig.defaultPythonPath = "python"

    val start = -5.0
    val stop = 5.0
    val pdfCount = 10000

    /*val estimatedPdf1D = DoubleKernelDensityEstimator.estimatePdf(
            samples = doubleArray1D(1000) { Random.nextGaussian() },
            start = start,
            stop = stop,
            pdfCount = pdfCount,
            kernel = ::gaussianKernel1D)*/

    /*showLine(title = "Kernel Density Estimation") {
        x = linearSpace(start, stop, pdfCount)
        y = estimatedPdf1D
    }*/

    val estimatedPdf2D = KernelDensityEstimator.estimatePdf(
            samples = doubleArray2D(2, 1000) { _, _ -> Random.nextGaussian() },
            start = D[-5.0, -5.0],
            stop = D[5.0, 5.0],
            pdfCount = 10000,
            kernel = ::gaussianKernel2D
    )

    showFigure {
        allPlots {
            position {
                rowCount = 1
                columnCount = 2
            }
        }

        plot {
            line {
                x = linearSpace(-5.0, 5.0, 10000)
                y = estimatedPdf2D.distribution[All,0]
            }
            position {
                column = 0
                row = 0
            }
        }

        plot {
            line {
                x = linearSpace(-5.0, 5.0, 10000)
                y = estimatedPdf2D.distribution[All,1]
            }
            position {
                column = 1
                row = 0
            }
        }
    }
}
