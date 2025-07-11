package repicea.math;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import repicea.stats.StatisticalUtility;

public abstract class AbstractComplexNumberSimulationStudy {

	/**
	 * A population unit.<p>
	 * The variable x follows a uniform distribution (3,10).
	 * @author Mathieu Fortin - December 2023
	 */
	static class PopulationUnit {
		private final double x;
		private final double y;
		@SuppressWarnings("unused")
		private final int id;
		private PopulationUnit(Matrix trueBeta, double trueStd, double x, int id) {
			this.id = id;
			this.x = x;
			y = trueBeta.getValueAt(0, 0) + x * trueBeta.getValueAt(1, 0) + StatisticalUtility.getRandom().nextGaussian() * trueStd;
		}
	}

	/**
	 * A sample of population units.
	 * @author Mathieu Fortin - December 2023
	 */
	@SuppressWarnings("serial")
	static class Sample extends ArrayList<PopulationUnit> {

		static final List<Double> X_VALUES = Arrays.asList(new Double[] {2d, 4d, 6d, 8d, 10d});

		static Sample createSample(Matrix trueBeta, double trueVariance, int sampleSize) {
			if (sampleSize%X_VALUES.size() != 0) {
				throw new InvalidParameterException("Requested sample size " + sampleSize + " is not a multiple of the number of x values: " + X_VALUES.size());
			}
			int nbRuns = sampleSize / X_VALUES.size();
			Sample s = new Sample();
			double trueStd = Math.sqrt(trueVariance);
			for (int i = 0; i < nbRuns; i++) {
				for (double x : X_VALUES) {
					s.add(new PopulationUnit(trueBeta, trueStd, x, i));
				}
			}
			return s;
		}

		
		private Matrix getMatrixX() {
			Matrix xMat = new Matrix(size(), 2);
			for (int i = 0; i < size(); i++) {
				xMat.setValueAt(i, 0, 1d);
				xMat.setValueAt(i, 1, get(i).x);
			}
			return xMat;
		}

		private Matrix getVectorY() {
			Matrix yVec = new Matrix(size(), 1);
			for (int i = 0; i < size(); i++) {
				yVec.setValueAt(i, 0, get(i).y);
			}
			return yVec;
		}
		
//		private DataSet convertIntoDataSet() {
//			DataSet ds = new DataSet(Arrays.asList(new String[] {"y", "x"}));
//			Object[] observation;
//			for (int i = 0; i < size(); i++) {
//				observation = new Object[2];
//				observation[0] = get(i).y;
//				observation[1] = get(i).x;
//				ds.addObservation(observation);
//			}
//			ds.indexFieldType();
//			return ds;
//		}
		
		
		Model getModel(AbstractComplexNumberSimulationStudy study) {
			Matrix xMat = getMatrixX();
			Matrix yVec = getVectorY();
			SymmetricMatrix invXtX = SymmetricMatrix.convertToSymmetricIfPossible(xMat.transpose().multiply(xMat).getInverseMatrix());
			Matrix betaHat = invXtX.multiply(xMat.transpose()).multiply(yVec);
			Matrix res = yVec.subtract(xMat.multiply(betaHat));
			int upsilon = yVec.m_iRows - xMat.m_iCols;
			double sigma2Hat = res.transpose().multiply(res).getValueAt(0, 0) / upsilon;
			return study.createModel(this, betaHat, invXtX, sigma2Hat, upsilon);
//			return new Model(this, betaHat, invXtX, sigma2Hat, upsilon);
		}
	}
	
	abstract Model createModel(Sample s, Matrix betaHat, SymmetricMatrix invXtX, double sigma2Hat, int upsilon);
	
	/**
	 * A model fitted to a sample of population units.
	 * @author Mathieu Fortin - December 2023
	 */
	static abstract class Model {
		final Matrix betaHat;
		final SymmetricMatrix invXtX;
		final double sigma2Hat;
		private Matrix invXtXChol;
		final int upsilon;
		
		Model(Sample sample, Matrix betaHat, SymmetricMatrix invXtX, double sigma2Hat, int upsilon) {
			this.betaHat = betaHat;
			this.invXtX = invXtX; 
			this.sigma2Hat = sigma2Hat;
			this.upsilon = upsilon;
		}
		
		Matrix getOmegaChol() {
			if (invXtXChol == null) {
				invXtXChol = invXtX.getLowerCholTriangle();
			}
			return invXtXChol;
		}
		
//		private Matrix getRandomDeviate(List<Double> xValues) {
//			Matrix xMat = createMatrixX(xValues);
//			double chiSquareDeviate = StatisticalUtility.getRandom().nextChiSquare(upsilon);
//			double sigma2_b = sigma2Hat + sigma2Hat * (chiSquareDeviate/upsilon - 1);
//			Matrix betaHat_b = betaHat.add(
//					getOmegaChol().multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN)).scalarMultiply(Math.sqrt(sigma2_b)));
//			Matrix result = xMat.multiply(betaHat_b).scalarAdd(0.5 * sigma2_b);
//			return result.expMatrix();
//		}
		
		static Matrix createMatrixX(List<Double> xValues) {
			Matrix xMat = new Matrix(xValues.size(),2);
			for (int i = 0; i < xValues.size(); i++) {
				xMat.setValueAt(i, 0, 1d);
				xMat.setValueAt(i, 1, xValues.get(i));
			}
			return xMat;
		}
		
		Matrix getXBeta(List<Double> xValues) {
			Matrix xMat = createMatrixX(xValues);
			return xMat.multiply(betaHat);
		}

//		private Matrix getBeauchampAndOlsonEstimator(List<Double> xValues) {
//			Matrix xMat = createMatrixX(xValues);
//			Matrix xBeta = getXBeta(xValues);
//			
//			Matrix output = new Matrix(xValues.size(), 1);
//			
//			for (int ii = 0; ii < xValues.size(); ii++) {
//				double term1 = Math.exp(xBeta.getValueAt(ii, 0) + sigma2Hat * .5);
//				Matrix xMat_ii = xMat.getSubMatrix(ii, ii, 0, xMat.m_iCols - 1);
//				int n = upsilon + xMat_ii.m_iCols;
//				Matrix xInvXtXxT = xMat_ii.multiply(invXtX).multiply(xMat_ii.transpose());
//				double phi = xInvXtXxT.getValueAt(0, 0) * n;
//				double sigma2hat2 = sigma2Hat * sigma2Hat;
//				double term2 = 1 - sigma2Hat * (2*phi + sigma2Hat)/(4*n) + 
//						sigma2hat2 * (sigma2hat2 + 2 * (16d/3 + 2 * phi) * sigma2Hat + 4 * phi * phi + 16 * phi) / (32 * n * n);
//				output.setValueAt(ii, 0, term1 * term2);
//			}
//			return output;
//		}

//		private Matrix getBaskervilleEstimator(List<Double> xValues) {
//			Matrix xBeta = getXBeta(xValues);
//			return xBeta.scalarAdd(sigma2Hat * .5).expMatrix();
//		}
//		
//		private ComplexMatrix getComplexRandomDeviate(List<Double> xValues) {
//			double chiSquareDeviate = StatisticalUtility.getRandom().nextChiSquare(upsilon);
//			ComplexNumber sigma2Hat_b = new ComplexNumber(sigma2Hat, sigma2Hat * (chiSquareDeviate/upsilon - 1));
//			ComplexNumber sigmaHat_b = sigma2Hat_b.sqrt();
//			Matrix xMat = createMatrixX(xValues);
//			Matrix betaDeviatesMat = getOmegaChol().multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN)); // here we get C * epsilon
//			
//			ComplexNumber[] cnArrayMean = new ComplexNumber[xValues.size()];
//			for (int ii = 0; ii < xValues.size(); ii++) {
//				Matrix xMat_ii = xMat.getSubMatrix(ii, ii, 0, xMat.m_iCols - 1);
//				ComplexNumber betaDeviates = new ComplexNumber(0, xMat_ii.multiply(betaDeviatesMat).getValueAt(0, 0)); // here we get 0 + x C epsilon i 
//				betaDeviates = betaDeviates.multiply(sigmaHat_b); // here we get 0 + sqrt{sigma2Hat_b} x C epsilon i
//				double xBetaHat = xMat_ii.multiply(betaHat).getValueAt(0, 0); 
//				ComplexNumber mean = betaDeviates.add(xBetaHat); // here we get x beta + + sqrt{sigma2Hat_b} x C epsilon i
//				ComplexNumber meanPlusCorrectionFactor = mean.add(sigma2Hat_b.multiply(0.5)); // here we get x beta + + sqrt{sigma2Hat_b} x C epsilon i + 0.5 * sigma2_b  				
//				cnArrayMean[ii] = meanPlusCorrectionFactor.exp();
//			}
//			return new ComplexMatrix(cnArrayMean); 
//		}
	}

}
