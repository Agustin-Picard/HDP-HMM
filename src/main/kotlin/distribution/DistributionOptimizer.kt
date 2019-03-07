package distribution

import tomasvolker.numeriko.core.interfaces.array1d.double.DoubleArray1D

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