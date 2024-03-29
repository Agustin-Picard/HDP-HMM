package distribution

import tomasvolker.kyplot.dsl.*
import tomasvolker.kyscript.KyScriptConfig
import tomasvolker.numeriko.core.dsl.D
import tomasvolker.numeriko.core.functions.*
import tomasvolker.numeriko.core.index.All
import tomasvolker.numeriko.core.interfaces.array1d.double.DoubleArray1D
import tomasvolker.numeriko.core.interfaces.array1d.double.elementWise
import tomasvolker.numeriko.core.interfaces.array2d.double.DoubleArray2D
import tomasvolker.numeriko.core.interfaces.factory.*
import tomasvolker.numeriko.core.linearalgebra.linearSpace
import tomasvolker.numeriko.core.operations.stack
import tomasvolker.numeriko.core.operations.unstack
import tomasvolker.numeriko.core.primitives.sumDouble
import utils.gaussianKernel1D
import utils.gaussianKernel2D
import utils.sumVector
import kotlin.math.pow
import kotlin.random.Random

typealias Kernel<T> = (T) -> T
typealias DoubleKernel = (Double) -> Double

interface DoubleDensityEstimator {
    fun estimatePdf(samples: Iterable<Number>,
                    start: Double,
                    stop: Double,
                    pdfCount: Int): SampledDistribution1D
}

interface DensityEstimator {
    fun estimatePdf(samples: DoubleArray2D,
                    start: DoubleArray1D,
                    stop: DoubleArray1D,
                    pdfCount: Int): SampledDistribution2D
}

class DoubleKernelDensityEstimator(val kernel: DoubleKernel): DoubleDensityEstimator {

    override fun estimatePdf(samples: Iterable<Number>,
                    start: Double,
                    stop: Double,
                    pdfCount: Int): SampledDistribution1D {
        val sampleArray = samples.toDoubleArray1D()
        val bandwidth = silvermanBandwidth(sampleArray)

        return SampledDistribution1D(
            distribution = doubleArray1D(pdfCount) { t ->
                sumDouble(0 until sampleArray.size) { i ->
                    kernel(((start + t * (stop - start) / pdfCount) - sampleArray[i]) / bandwidth)
                }
            } / (bandwidth * sampleArray.size),
            start = start,
            stop = stop
        )
    }

    private fun silvermanBandwidth(samples: DoubleArray1D) =
        (4.0 * samples.std() / (3.0 * samples.size)).pow(0.2)

}


class KernelDensityEstimator(val kernel: Kernel<DoubleArray1D>): DensityEstimator {

    override fun estimatePdf(
        samples: DoubleArray2D,
        start: DoubleArray1D,
        stop: DoubleArray1D,
        pdfCount: Int
    ): SampledDistribution2D {
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


fun main() {
    KyScriptConfig.defaultPythonPath = "python"

    val start = -5.0
    val stop = 5.0
    val pdfCount = 10000

    val estimatedPdf1D = DoubleKernelDensityEstimator(::gaussianKernel1D).estimatePdf(
            samples = doubleArray1D(1000) { Random.nextGaussian() },
            start = start,
            stop = stop,
            pdfCount = pdfCount)

    showLine(title = "Kernel Density Estimation") {
        x = linearSpace(start, stop, pdfCount)
        y = estimatedPdf1D.distribution
    }

    /*val estimatedPdf2D = KernelDensityEstimator(::gaussianKernel2D).estimatePdf(
            samples = doubleArray2D(2, 1000) { _, _ -> Random.nextGaussian() },
            start = D[-5.0, -5.0],
            stop = D[5.0, 5.0],
            pdfCount = 10000
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
    }*/
}
