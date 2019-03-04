package hdphmm

import tomasvolker.numeriko.core.interfaces.array1d.double.DoubleArray1D
import tomasvolker.numeriko.core.interfaces.array2d.double.DoubleArray2D
import tomasvolker.numeriko.core.primitives.sumDouble
import kotlin.math.ln


typealias SampledDistribution1D = SampledDistribution<Double,DoubleArray1D>
typealias SampledDistribution2D = SampledDistribution<DoubleArray1D,DoubleArray2D>

data class SampledDistribution<P,T>(val distribution: T, val start: P, val stop: P)

inline val <reified P> SampledDistribution<P,*>.sampleCount get() =
    if (distribution is DoubleArray2D) distribution.shape1
    else if (distribution is DoubleArray1D) distribution.size
    else throw IllegalArgumentException("Invalid distribution format")

fun DoubleDistribution<*>.estimateKLDivergence(other: SampledDistribution1D): Double =
        pdf(nSamples = other.sampleCount,
                start = other.start,
                stop = other.stop).let {
            sumDouble(0 until other.sampleCount) { i -> it[i] * ln(other.distribution[i] / it[i]) }
        }