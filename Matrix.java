package MatrixMath;

import java.util.Arrays;

public class Matrix {
	private double[][] matrix;
	private int rows, columns;

	public Matrix(double[][] matrix) throws InvalidMatrixException {
		if(matrix.length != matrix[0].length) {
			throw new InvalidMatrixException();
		}
		
		for(int i=0; i<matrix.length-1; ++i) {
			if(matrix[i].length != matrix[i+1].length) {
				throw new InvalidMatrixException();
			}
		}
		
		this.matrix = matrix;
		this.rows = matrix.length;
		this.columns = matrix[0].length;
	}
	
	public static Matrix add(Matrix m1, Matrix m2) throws MatrixAdditiveException, InvalidMatrixException {
		if(m1.rows != m2.rows || m1.columns != m2.columns) {
			throw new MatrixAdditiveException();
		}
		
		double[][] ans = new double[m1.rows][m1.columns];
		
		for(int i=0; i<m1.rows; ++i) {
			for(int j=0; j<m1.columns; ++j) {
				ans[i][j] = m1.matrix[i][j] + m2.matrix[i][j];
			}
		}
		
		return new Matrix(ans);
	}
	public static Matrix subtract(Matrix m1, Matrix m2) throws MatrixAdditiveException, InvalidMatrixException {
		if(m1.rows != m2.rows || m1.columns != m2.columns) {
			throw new MatrixAdditiveException();
		}
		
		double[][] ans = new double[m1.rows][m1.columns];
		
		for(int i=0; i<m1.rows; ++i) {
			for(int j=0; j<m1.columns; ++j) {
				ans[i][j] = m1.matrix[i][j] - m2.matrix[i][j];
			}
		}
		
		return new Matrix(ans);
	}
	public static Matrix multiplyM(Matrix m1, Matrix m2) throws MatrixMultiplicativeException, InvalidMatrixException {
		if(m1.columns != m2.rows) {
			throw new MatrixMultiplicativeException();
		}
		
		double[][] ans = new double[m1.rows][m2.columns];
		
		for (int i=0; i<m1.rows; ++i) {
            for (int j=0; j < m2.columns; ++j) {
                for (int k=0; k<m1.columns; ++k) {
                    ans[i][j] += m1.matrix[i][k] * m2.matrix[k][j];
                }
            }
        }
		
		return new Matrix(ans);
	}
	public static Matrix divideM(Matrix m1, Matrix m2) throws MatrixMultiplicativeException, InvalidMatrixException, SquareMatrixException, MatrixComputationTimeException, InvertibleMatrixException {
		Matrix invMatrix2 = inverse(m2);
		return multiplyM(m1, invMatrix2);
	}
	public static Matrix inverse(Matrix m1) throws SquareMatrixException, MatrixComputationTimeException, InvalidMatrixException, InvertibleMatrixException {
		if(m1.rows != m1.columns) {
			throw new SquareMatrixException();
		}
		
		double det = getDet(m1);
		if(det==0) {
			throw new InvertibleMatrixException();
		}
		Matrix adj = adjoint(m1);
		
		return multiplyS(adj, det);
	}
	public static Matrix adjoint(Matrix m1) throws SquareMatrixException, MatrixComputationTimeException, InvalidMatrixException {
		if(m1.rows != m1.columns) {
			throw new SquareMatrixException();
		}
		if(m1.rows!=2) {
			throw new MatrixComputationTimeException();
		}
		
		double[][] ans = new double[2][2];
		
		ans[0][0] = m1.matrix[1][1];
		ans[1][1] = m1.matrix[0][0];
		ans[0][1] = -m1.matrix[0][1];
		ans[1][0] = -m1.matrix[1][0];
		
		return new Matrix(ans);
	}
	public static Matrix transpose(Matrix m1) throws InvalidMatrixException {
		double[][] ans = new double[m1.columns][m1.rows];
		for(int i=0; i<m1.rows; ++i) {
			for(int j=0; j<m1.columns; ++j) {
				ans[j][i] = m1.matrix[i][j];
			}
		}
		return new Matrix(ans);
	}
	public static Matrix multiplyS(Matrix m1, double scalar) throws InvalidMatrixException {
		double[][] ans = new double[m1.rows][m1.columns];
		
		for(int i=0; i<m1.rows; ++i) {
			for(int j=0; j<m1.columns; ++j) {
				ans[i][j] = scalar * m1.matrix[i][j];
			}
		}
		
		return new Matrix(ans);
	}
	public static Matrix divideS(Matrix m1, double scalar) throws InvalidMatrixException {
		double[][] ans = new double[m1.rows][m1.columns];
		
		for(int i=0; i<m1.rows; ++i) {
			for(int j=0; j<m1.columns; ++j) {
				ans[i][j] = 1/scalar * m1.matrix[i][j];
			}
		}
		
		return new Matrix(ans);
	}
	
	public static double getDet(Matrix m1) throws SquareMatrixException, MatrixComputationTimeException {
		if(m1.rows != m1.columns) {
			throw new SquareMatrixException();
		}
		if(m1.rows!=2) {
			throw new MatrixComputationTimeException();
		}
		if(m1.rows==1) return m1.matrix[0][0];
		
				
		return m1.matrix[0][0]*m1.matrix[1][1] - m1.matrix[0][1]*m1.matrix[1][0];
	}
	
	@Override
	public String toString() {
		return Arrays.deepToString(this.matrix);
	}
}





class InvalidMatrixException extends Exception {
	public InvalidMatrixException() {
		super();
		System.out.println("This matrix is not valid for this operation.");
	}
}

class InvertibleMatrixException extends Exception {
	public InvertibleMatrixException() {
		super();
		System.out.println("This matrix is not invertible because it has a determinant of zero.");
	}
}

class MatrixAdditiveException extends Exception {
	public MatrixAdditiveException() {
		super();
		System.out.println("These matrices cannot be added together due to differences in their order/dimensions.");
	}
}

class MatrixMultiplicativeException extends Exception {
	public MatrixMultiplicativeException() {
		super();
		System.out.println("Multiplication between these two matrices cannot be performed due to inequality in the matrices' rows and columns.");
	}
}

class SquareMatrixException extends Exception {
	public SquareMatrixException() {
		super();
		System.out.println("This operation can only be perfomed on square matrices (a matrix where the number of rows is equal to the number of columns).");
	}
}

class MatrixComputationTimeException extends Exception {
	public MatrixComputationTimeException() {
		super();
		System.out.println("Performing this operation takes a LOT of computational time and power. As such, this program does not support performing this operation on matrices of this size.");
	}
}
