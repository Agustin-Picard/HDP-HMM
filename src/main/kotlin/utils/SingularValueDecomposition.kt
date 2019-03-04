package utils

import koma.extensions.fill
import koma.extensions.get
import koma.matrix.Matrix
import koma.zeros
import tomasvolker.numeriko.core.interfaces.array2d.double.DoubleArray2D
import tomasvolker.numeriko.core.interfaces.factory.doubleArray2D

fun Matrix<Double>.toDoubleArray2D() =
    doubleArray2D(numRows(), numCols()) { i, j -> this[i, j] }

fun DoubleArray2D.toMatrixDouble() = zeros(shape0, shape1).fill { i, j -> this[i, j] }

data class SingularValueDecomposition(
    val u: DoubleArray2D,
    val s: DoubleArray2D,
    val v: DoubleArray2D
)

fun DoubleArray2D.singularValueDecomposition(): SingularValueDecomposition {
    val (u, s, v) = toMatrixDouble().SVD()
    return SingularValueDecomposition(
        u = u.toDoubleArray2D(),
        s = s.toDoubleArray2D(),
        v = v.toDoubleArray2D()
    )
}
