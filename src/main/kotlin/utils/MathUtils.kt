package utils

import hdphmm.mean
import org.apache.commons.math3.special.Erf
import org.apache.commons.math3.special.Gamma
import tomasvolker.numeriko.core.functions.average
import tomasvolker.numeriko.core.functions.exp
import tomasvolker.numeriko.core.functions.matMul
import tomasvolker.numeriko.core.functions.transpose
import tomasvolker.numeriko.core.interfaces.array1d.double.DoubleArray1D
import tomasvolker.numeriko.core.interfaces.array1d.double.elementWise
import tomasvolker.numeriko.core.interfaces.array2d.double.DoubleArray2D
import tomasvolker.numeriko.core.interfaces.array2d.double.elementWise
import tomasvolker.numeriko.core.primitives.squared
import tomasvolker.numeriko.core.primitives.sumDouble
import kotlin.math.PI
import kotlin.math.exp
import kotlin.math.pow
import kotlin.math.sqrt

fun gammaFunction(x: Double): Double =
        Gamma.gamma(x)

fun logGammaFunction(x: Double) =
        Gamma.logGamma(x)

fun erfinvFunction(x: Double) =
        Erf.erfInv(x)

fun Int.factorial(): Int =
    if (this > 0) this * factorial() else 1

fun sphericalKernel1D(x: Double) = if (x <= 1.0) 1.0 else 0.0

fun DoubleArray1D.estimateKurtosis() =
        mean().let { mean ->
            sumDouble(0 until size) { i -> (this[i] - mean).pow(4) } /
                    (sumDouble(0 until size) { i -> (this[i] - mean).squared() }.squared() / size) - 3.0
        }

fun DoubleArray2D.pow(x: Double) =
        singularValueDecomposition().let {
            it.u matMul it.s.elementWise { it.pow(x) } matMul it.v.transpose()
        }

fun DoubleArray1D.unbiased() = this - average()

fun gaussianKernel1D(x: Double) = exp(-(x * x) / 2.0) / sqrt(2.0 * PI)

fun gaussianKernel2D(x: DoubleArray1D) =
        exp(x.elementWise { -it.squared() / 2.0 }) / sqrt(
            (2.0 * PI).pow(
                x.size
            )
        )

fun sumVector(indices: IntProgression, selector: (Int) -> DoubleArray1D) =
        indices.asSequence()
                .map { selector(it) }
                .reduce { acc, next -> acc + next }