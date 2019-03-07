package utils

import distribution.ExponentialDistribution
import distribution.GammaDistribution
import distribution.GaussianDistribution1D
import distribution.LogNormalDistribution
import kotlin.random.Random

object TraceGenerator {

    fun gaussianTrace(mean: Double, std: Double, nSamples: Int, random: Random = Random.Default) =
            GaussianDistribution1D(mean, std).nextSample(nSamples, random)

    fun exponentialTrace(rate: Double, nSamples: Int, random: Random = Random.Default) =
            ExponentialDistribution(rate).nextSample(nSamples, random)

    fun logNormalTrace(mean: Double, std: Double, nSamples: Int, random: Random = Random.Default) =
            LogNormalDistribution(mean, std).nextSample(nSamples, random)

    fun gammaTrace(shape: Double, scale: Double, nSamples: Int, random: Random = Random.Default) =
            GammaDistribution(shape, scale).nextSample(nSamples, random)
}