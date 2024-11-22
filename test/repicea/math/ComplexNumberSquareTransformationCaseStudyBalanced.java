/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2023 His Majesty the King in Right of Canada
 * Author: Mathieu Fortin, Canadian Wood Fibre Centre, Canadian Forest Service
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 *
 * This library is distributed with the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * Please see the license at http://www.gnu.org/copyleft/lesser.html.
 */
package repicea.math;

import java.io.File;
import java.io.IOException;
import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import repicea.io.FormatField;
import repicea.io.javacsv.CSVField;
import repicea.io.javacsv.CSVWriter;
import repicea.lang.REpiceaSystem;
import repicea.stats.Distribution;
import repicea.stats.StatisticalUtility;
import repicea.stats.estimates.ComplexMonteCarloEstimate;
import repicea.stats.estimates.ConfidenceInterval;
import repicea.stats.estimates.MonteCarloEstimate;

/**
 * The ComplexNumberSimpleCaseStudy class implements a simulation study around 
 * the square transformation of Gaussian variables.
 * @author Mathieu Fortin - Oct 2023
 */
public class ComplexNumberSquareTransformationCaseStudyBalanced {

	/**
	 * A population unit.<p>
	 * The variable x follows a uniform distribution (3,10).
	 * @author Mathieu Fortin - December 2023
	 */
	private static class PopulationUnit {
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
	private static class Sample extends ArrayList<PopulationUnit> {

		static final List<Double> X_VALUES = Arrays.asList(new Double[] {1d, 2d, 3d, 4d, 5d, 6d, 7d, 8d, 9d, 10d});
		
		private static Sample createSample(Matrix trueBeta, double trueVariance, int sampleSize) {
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
		
		
		private Model getModel() {
			Matrix xMat = getMatrixX();
			Matrix yVec = getVectorY();
			SymmetricMatrix invXtX = SymmetricMatrix.convertToSymmetricIfPossible(xMat.transpose().multiply(xMat).getInverseMatrix());
			Matrix betaHat = invXtX.multiply(xMat.transpose()).multiply(yVec);
			Matrix res = yVec.subtract(xMat.multiply(betaHat));
			int upsilon = yVec.m_iRows - xMat.m_iCols;
			double sigma2Hat = res.transpose().multiply(res).getValueAt(0, 0) / upsilon;
			return new Model(this, betaHat, invXtX, sigma2Hat, upsilon);
		}
		
	}
	
	/**
	 * A model fitted to a sample of population units.
	 * @author Mathieu Fortin - December 2023
	 */
	private static class Model {
		private final Matrix betaHat;
		private final SymmetricMatrix invXtX;
		private final double sigma2Hat;
		private Matrix invXtXChol;
		private final int upsilon;
		
		Model(Sample sample, Matrix betaHat, SymmetricMatrix invXtX, double sigma2Hat, int upsilon) {
			this.betaHat = betaHat;
			this.invXtX = invXtX; 
			this.sigma2Hat = sigma2Hat;
			this.upsilon = upsilon;
		}
		
		private Matrix getOmegaChol() {
			if (invXtXChol == null) {
				invXtXChol = invXtX.getLowerCholTriangle();
			}
			return invXtXChol;
		}
		
		private Matrix getRandomDeviate(List<Double> xValues) {
			Matrix xMat = createMatrixX(xValues);
			double chiSquareDeviate = StatisticalUtility.getRandom().nextChiSquare(upsilon);
			double sigma2_b = sigma2Hat + sigma2Hat * (chiSquareDeviate/upsilon - 1);
			Matrix betaHat_b = betaHat.add(
					getOmegaChol().multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN)).scalarMultiply(Math.sqrt(sigma2_b)));
			Matrix result = xMat.multiply(betaHat_b).elementWisePower(2d).scalarAdd(sigma2_b);
			return result;
		}
		
		private static Matrix createMatrixX(List<Double> xValues) {
			Matrix xMat = new Matrix(xValues.size(),2);
			for (int i = 0; i < xValues.size(); i++) {
				xMat.setValueAt(i, 0, 1d);
				xMat.setValueAt(i, 1, xValues.get(i));
			}
			return xMat;
		}
		
//		private Matrix getXBeta(List<Double> xValues) {
//			Matrix xMat = createMatrixX(xValues);
//			return xMat.multiply(betaHat);
//		}

		private Matrix getGregoireEstimator(List<Double> xValues) {
			Matrix xMat = createMatrixX(xValues);
			Matrix varMuHat = xMat.multiply(invXtX).multiply(xMat.transpose()).scalarMultiply(sigma2Hat);
			Matrix xBeta = xMat.multiply(betaHat);
			Matrix pred = xBeta.elementWisePower(2d).scalarAdd(sigma2Hat).add(varMuHat.diagonalVector().scalarMultiply(-1d));
			return pred;
		}

		private Matrix getGregoireVarianceEstimator(List<Double> xValues) {
			Matrix xMat = createMatrixX(xValues);
			Matrix varMuHat = xMat.multiply(invXtX).multiply(xMat.transpose()).scalarMultiply(sigma2Hat);
			Matrix xBeta = xMat.multiply(betaHat);
			Matrix varianceHat = xBeta.elementWisePower(2d).elementWiseMultiply(varMuHat.diagonalVector()).scalarMultiply(4d).
					subtract(varMuHat.diagonalVector().elementWisePower(2d)).scalarAdd(2*sigma2Hat*sigma2Hat/upsilon);
			return varianceHat;
		}

		private ComplexMatrix getComplexRandomDeviate(List<Double> xValues) {
			double chiSquareDeviate = StatisticalUtility.getRandom().nextChiSquare(upsilon);
			ComplexNumber sigma2Hat_b = new ComplexNumber(sigma2Hat, sigma2Hat * (chiSquareDeviate/upsilon - 1));
			ComplexNumber sigmaHat_b = sigma2Hat_b.sqrt();
			Matrix xMat = createMatrixX(xValues);
			Matrix betaDeviatesMat = getOmegaChol().multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN)); // here we get C * epsilon
			
			ComplexNumber[] cnArrayMean = new ComplexNumber[xValues.size()];
			for (int ii = 0; ii < xValues.size(); ii++) {
				Matrix xMat_ii = xMat.getSubMatrix(ii, ii, 0, xMat.m_iCols - 1);
				ComplexNumber betaDeviates = new ComplexNumber(0, xMat_ii.multiply(betaDeviatesMat).getValueAt(0, 0)); // here we get 0 + x C epsilon i 
				betaDeviates = betaDeviates.multiply(sigmaHat_b); // here we get 0 + sqrt{sigma2Hat_b} x C epsilon i
				double xBetaHat = xMat_ii.multiply(betaHat).getValueAt(0, 0); 
				ComplexNumber mean = betaDeviates.add(xBetaHat); // here we get x beta + + sqrt{sigma2Hat_b} x C epsilon i
				ComplexNumber meanPlusCorrectionFactor = mean.pow(2).add(sigma2Hat_b); // here we get (x beta + sqrt{sigma2Hat_b} x C epsilon i)^2 + sigma2_b  				
				cnArrayMean[ii] = meanPlusCorrectionFactor;
			}
			return new ComplexMatrix(cnArrayMean); 
		}
	}
	
//	private static CSVWriter createWriterForIndividualEstimate(String filename) throws IOException {
//		CSVWriter writer = new CSVWriter(new File(filename.concat("exampleSingleEstimate.csv")), false);
//		List<FormatField> fields = new ArrayList<FormatField>();
//		fields.add(new CSVField("RealID"));
//		fields.add(new CSVField("x"));
//		fields.add(new CSVField("real"));
//		fields.add(new CSVField("imag"));
//		writer.setFields(fields);
//		return writer;
//	}
	
	private void doRun(int sampleSize, 
			int nbRealizations, 
			double b0, 
			double b1, 
			double s2,
			String filename) throws IOException {
		System.out.println("Simulating [" + b0 + "; " + b1 + "; " + s2 +"] with sample size n = " + sampleSize);
		int nbInnerReal = 10000;
		CSVWriter writer = new CSVWriter(new File(filename.concat("" + sampleSize + ".csv")), false);
		List<FormatField> fields = new ArrayList<FormatField>();
		fields.add(new CSVField("RealID"));
		fields.add(new CSVField("b0"));
		fields.add(new CSVField("b1"));
		fields.add(new CSVField("sigma2"));
		fields.add(new CSVField("x"));
		fields.add(new CSVField("MC_meanEst"));
		fields.add(new CSVField("CMC_meanEst"));
		fields.add(new CSVField("CMC_varEst"));
		fields.add(new CSVField("Lower95"));
		fields.add(new CSVField("Upper95"));
		fields.add(new CSVField("Lower99"));
		fields.add(new CSVField("Upper99"));
		fields.add(new CSVField("GregoireEstimator"));
		fields.add(new CSVField("GregoireVarianceEstimator"));
		fields.add(new CSVField("TrueExp"));
		fields.add(new CSVField("TrueMSE"));
		writer.setFields(fields);
		
		Matrix trueBeta = new Matrix(2,1);
		trueBeta.setValueAt(0, 0, b0);
		trueBeta.setValueAt(1, 0, b1);
		double trueVariance = s2;
		
		
		List<Double> xValues = Sample.X_VALUES;
		Matrix mse = null;
		Matrix theta = null;
		for (int real = 1; real <= nbRealizations; real++) {
			if (real % 1000 == 0) {
				System.out.println("Running realization " + real);
			}

			Sample s = Sample.createSample(trueBeta, trueVariance, sampleSize);
			Model m = s.getModel();
						
			if (real == 1) {
				mse  = new Matrix(xValues.size(),1);
				theta  = new Matrix(xValues.size(),1);
				Matrix xMat = Model.createMatrixX(xValues);
				Matrix Omega = xMat.multiply(m.invXtX).multiply(xMat.transpose());
				Matrix xBeta = xMat.multiply(trueBeta);
				for (int i = 0; i < xValues.size(); i++) {
					double xBeta2 = xBeta.getValueAt(i, 0) * xBeta.getValueAt(i, 0);
					double omega = Omega.getValueAt(i, i);
					double mseValue = 4 * xBeta2 * omega * trueVariance + 
							2 * trueVariance * trueVariance * (omega * omega + (1-omega)*(1-omega)/m.upsilon);
					mse.setValueAt(i, 0, mseValue);
					theta.setValueAt(i, 0, xBeta2 + trueVariance);
				}
			}
			
			MonteCarloEstimate mcEstimator = new MonteCarloEstimate();
			ComplexMonteCarloEstimate cmcEstimator = new ComplexMonteCarloEstimate();
//			CSVWriter singleEstWriter = null;
//			if (real == 1) {
//				singleEstWriter = createWriterForIndividualEstimate(filename);
//			}
			
			for (int innerReal = 0; innerReal < nbInnerReal; innerReal++) {
				mcEstimator.addRealization(m.getRandomDeviate(xValues)); 
				ComplexMatrix complexRealizations = m.getComplexRandomDeviate(xValues);
				cmcEstimator.addRealization(complexRealizations); 
//				if (real == 1) {
//					for (int innerLoc = 0; innerLoc < complexRealizations.m_iRows; innerLoc++) {
//						Object[] record = new Object[4];
//						record[0] = innerReal;
//						record[1] = xValues.get(innerLoc);
//						record[2] = complexRealizations.getValueAt(innerLoc, 0).realPart;
//						record[3] = complexRealizations.getValueAt(innerLoc, 0).imaginaryPart;
//						singleEstWriter.addRecord(record);
//					}
//				}
			}
//			if (real == 1) {
//				singleEstWriter.close();
//			}
			
			Matrix mcMean = mcEstimator.getMean();
			ComplexMatrix cmcMean = cmcEstimator.getMean();
			ComplexSymmetricMatrix cmcPsVar = cmcEstimator.getPseudoVariance();
			
			ComplexMonteCarloEstimate imaginaryRealizations = cmcEstimator.getImaginaryRealizations();
			ConfidenceInterval ci95 = imaginaryRealizations.getConfidenceIntervalBounds(0.95);
			ConfidenceInterval ci99 = imaginaryRealizations.getConfidenceIntervalBounds(0.99);

			Matrix gregoireEstimator = m.getGregoireEstimator(xValues);
			Matrix gregoireVarianceEstimator = m.getGregoireVarianceEstimator(xValues);
			
			for (int ii = 0; ii < xValues.size(); ii++) {
				Object[] record = new Object[fields.size()];
				record[0] = real;
				record[1] = m.betaHat.getValueAt(0, 0);
				record[2] = m.betaHat.getValueAt(1, 0);
				record[3] = m.sigma2Hat;
				record[4] = xValues.get(ii);
				record[5] = mcMean.getValueAt(ii, 0);
				record[6] = cmcMean.getValueAt(ii, 0).realPart;
				record[7] = -cmcPsVar.getValueAt(ii, ii).realPart;
				record[8] = ci95.getLowerLimit().getValueAt(ii, 0);
				record[9] = ci95.getUpperLimit().getValueAt(ii, 0);
				record[10] = ci99.getLowerLimit().getValueAt(ii, 0);
				record[11] = ci99.getUpperLimit().getValueAt(ii, 0);
				record[12] = gregoireEstimator.getValueAt(ii, 0);
				record[13] = gregoireVarianceEstimator.getValueAt(ii, 0);
				record[14] = theta.getValueAt(ii, 0);
				record[15] = mse.getValueAt(ii, 0);
				writer.addRecord(record);
			}
		}
		writer.close();
	}
	
	public static void main(String[] args) throws IOException, InterruptedException {
		String filename = REpiceaSystem.retrieveArgument("-outdir", Arrays.asList(args));
		System.out.println("Export filename: " + filename);
		List<Double> variances = new ArrayList<Double>();
//		variances.add(0.5);
//		variances.add(1.0);
//		variances.add(2.0);
		variances.add(4.0);
		for (double variance : variances) {
			String rootFilename = filename.concat("caseStudySquareBalanced_" + variance + "_");
			List<Integer> sampleSizes = new ArrayList<Integer>();
			sampleSizes.add(20);
			sampleSizes.add(50);
			sampleSizes.add(100);
			sampleSizes.add(200);
			List<Thread> threads = new ArrayList<Thread>();
			for (int sampleSize : sampleSizes) {
				Runnable doRun = new Runnable() {
					public void run() {
						ComplexNumberSquareTransformationCaseStudyBalanced s = new ComplexNumberSquareTransformationCaseStudyBalanced();
						try {
							s.doRun(sampleSize, 50000, 5, 0.25, variance, rootFilename);
						} catch (IOException e) {
							e.printStackTrace();
						}
					};
				};
				Thread t = new Thread(doRun);
				threads.add(t);
				t.start();
			}
			for (Thread t : threads) {
				t.join();
			}
		}
	}	
}
