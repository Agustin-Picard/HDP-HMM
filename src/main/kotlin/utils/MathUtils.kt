package utils

import org.apache.commons.math3.special.Erf
import org.apache.commons.math3.special.Gamma
import tomasvolker.numeriko.core.functions.exp
import tomasvolker.numeriko.core.interfaces.array1d.double.DoubleArray1D
import tomasvolker.numeriko.core.interfaces.array1d.double.elementWise
import tomasvolker.numeriko.core.primitives.squared
import kotlin.math.PI
import kotlin.math.exp
import kotlin.math.pow
import kotlin.math.sqrt

fun gammaFunction(x: Double): Double =
        Gamma.gamma(x)

fun logGammaFunction(x: Double) =
        Gamma.logGamma(x)

fun digammaFunction(x: Double) =
        Gamma.digamma(x)

fun trigammaFunction(x: Double) =
        Gamma.trigamma(x)

fun erfinvFunction(x: Double) =
        Erf.erfInv(x)

fun Int.factorial(): Int =
    if (this > 0) this * factorial() else 1

fun sphericalKernel1D(x: Double) = if (x <= 1.0) 1.0 else 0.0

fun gaussianKernel1D(x: Double) = exp(-(x * x) / 2.0) / sqrt(2.0 * PI)

fun gaussianKernel2D(x: DoubleArray1D) =
        exp(x.elementWise { -it.squared() / 2.0 }) / sqrt(
            (2.0 * PI).pow(
                x.size
            )
        )

