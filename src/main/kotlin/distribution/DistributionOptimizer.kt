package distribution

import org.openrndr.color.ColorRGBa
import org.openrndr.math.Vector2
import tomasvolker.kyscript.KyScriptConfig
import tomasvolker.numeriko.core.interfaces.array1d.double.DoubleArray1D
import tomasvolker.numeriko.core.interfaces.factory.toDoubleArray1D
import tomasvolker.openrndr.math.plot.plotLine
import tomasvolker.openrndr.math.plot.quickPlot2D
import tomasvolker.openrndr.math.primitives.d
import utils.SegmentedTracesReader
import utils.gaussianKernel1D
import utils.unbiased

object DistributionOptimizer {

    fun optimize(data: DoubleArray1D, kernel: DoubleKernel, resolution: Int = 10000): DoubleDistribution<*> {
        val dataMin = data.min() ?: error("empty data array")
        val dataMax = data.max() ?: error("empty data array")
        val margin = (dataMax - dataMin) / 10.0
        val start = dataMin - margin
        val stop = dataMax + margin

        // calculate the empiric density through kde
        val empiricDensity = DoubleKernelDensityEstimator(kernel).estimatePdf(data, start, stop, resolution)

        // fit the parameters of the distributions we want to test
        val gaussianParameters = GaussianDistributionFitter1D.fitParameters(data)
        val logNormalParameters = LogNormalDistributionFitter.fitParameters(data)
        val exponentialParameters = ExponentialDistributionFitter.fitParameters(data)
        val gammaParameters = GammaDistributionFitter.fitParameters(data)

        // generate the corresponding distributions
        val normalDistribution = GaussianDistribution1D.fromParameters(gaussianParameters)
        val logNormalDistribution = LogNormalDistribution.fromParameters(logNormalParameters)
        val exponentialDistribution = ExponentialDistribution.fromParameters(exponentialParameters)
        val gammaDistribution = GammaDistribution.fromParameters(gammaParameters)

        // compute the kullback-leibler convergence
        val normalKLD = normalDistribution.estimateKLDivergence(empiricDensity)
        val logNormalKLD = logNormalDistribution.estimateKLDivergence(empiricDensity)
        val exponentialKLD = exponentialDistribution.estimateKLDivergence(empiricDensity)
        val gammaKLD = gammaDistribution.estimateKLDivergence(empiricDensity)

        return listOf<Pair<Double,DoubleDistribution<*>>>(
            normalKLD to normalDistribution,
            logNormalKLD to logNormalDistribution,
            exponentialKLD to exponentialDistribution,
            gammaKLD to gammaDistribution
        ).maxBy { it.first }?.second ?: error("Error in distributions")
    }
}

fun main() {
    KyScriptConfig.defaultPythonPath = "python"

    val filename = "data/msm_18725407_746.csv"
    val traces = SegmentedTracesReader(filename).readTraces()
    val segmentedTraces = SegmentedTracesReader(filename).readSegmented()
    val densityList = mutableListOf<SampledDistribution1D>()
    val distributionList = mutableListOf<DoubleDistribution<*>>()

    segmentedTraces.forEach {
        val data = it.data.toDoubleArray1D()
        val dataMin = data.min() ?: error("empty data array")
        val dataMax = data.max() ?: error("empty data array")
        val margin = (dataMax - dataMin) / 10.0
        val start = dataMin - margin
        val stop = dataMax + margin

        densityList.add(DoubleKernelDensityEstimator(::gaussianKernel1D).estimatePdf(data, start, stop, 1000))
        distributionList.add(DistributionOptimizer.optimize(data, ::gaussianKernel1D, 1000))
    }

    println(distributionList)
    println(distributionList[0].parameters)

    val start = densityList[0].start.also { println(it) }
    val stop = densityList[0].stop.also { println(it) }

    val plotDist = distributionList[0].pdf(1000, 1.0, 1000.0).also { println(it) }
        .mapIndexed { index, d -> Vector2(index.d / 10.0, d * 100.0) }
    val plotDensity = densityList[1].distribution.mapIndexed { index, d -> Vector2(index.d / 10.0, d * 100.0) }

    quickPlot2D {
        stroke = ColorRGBa.RED
        plotLine(plotDensity)
//        stroke = ColorRGBa.BLUE
//        plotLine(plotDist)
    }

    /*showFigure {
        allPlots {
            position {
                rowCount = 1
                columnCount = 2
            }
        }

        plot {
            line {
                x = linearSpace(plotStart, plotStop, 1000)
                y = densityList[1].distribution
            }
            position {
                column = 0
                row = 0
            }
        }

        plot {
            line {
                x = linearSpace(plotStart, plotStop, 1000)
                y = plotDist
            }
            position {
                column = 1
                row = 0
            }
        }
    }*/
}