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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import repicea.io.FormatField;
import repicea.io.javacsv.CSVField;
import repicea.io.javacsv.CSVWriter;
import repicea.math.utility.GammaUtility;
import repicea.math.utility.MathUtility;
import repicea.stats.Distribution;
import repicea.stats.Distribution.Type;
import repicea.stats.StatisticalUtility;
import repicea.stats.estimates.ComplexMonteCarloEstimate;
import repicea.stats.estimates.MonteCarloEstimate;
import repicea.stats.sampling.SamplingUtility;
import repicea.util.ObjectUtility;

/**
 * The ComplexNumberSimpleCaseStudy class implements a simulation study around 
 * the exponential, square and log transformation of Gaussian variables.
 * @author Mathieu Fortin - Oct 2023
 */
public class ComplexNumberSimpleCaseStudy {

	
	static class PopulationUnit {
		final double x;
		final double y;
		final int id;
		PopulationUnit(Matrix trueBeta, double trueStd, int id) {
			this.id = id;
			x = StatisticalUtility.getRandom().nextDouble() * 7 + 3;
			y = trueBeta.getValueAt(0, 0) + x * trueBeta.getValueAt(1, 0) + StatisticalUtility.getRandom().nextGaussian() * trueStd;
		}
	}

	
	@SuppressWarnings("serial")
	static class Sample extends ArrayList<PopulationUnit> {
		
		static Sample createSample(Matrix trueBeta, double trueVariance, int sampleSize) {
			Sample s = new Sample();
			double trueStd = Math.sqrt(trueVariance);
			for (int i = 0; i < sampleSize; i++) 
				s.add(new PopulationUnit(trueBeta, trueStd, i));
			return s;
		}
		
		Matrix getMatrixX() {
			Matrix xMat = new Matrix(size(), 2);
			for (int i = 0; i < size(); i++) {
				xMat.setValueAt(i, 0, 1d);
				xMat.setValueAt(i, 1, get(i).x);
			}
			return xMat;
		}

		Matrix getVectorY() {
			Matrix yVec = new Matrix(size(), 1);
			for (int i = 0; i < size(); i++) {
				yVec.setValueAt(i, 0, get(i).y);
			}
			return yVec;
		}
		
		Model getModel(Double sigma2) {
			Matrix xMat = getMatrixX();
			Matrix yVec = getVectorY();
			SymmetricMatrix invXtX = SymmetricMatrix.convertToSymmetricIfPossible(xMat.transpose().multiply(xMat).getInverseMatrix());
			Matrix betaHat = invXtX.multiply(xMat.transpose()).multiply(yVec);
			Matrix res = yVec.subtract(xMat.multiply(betaHat));
			int upsilon = yVec.m_iRows - xMat.m_iCols;
			double sigma2Hat = sigma2 != null ?
					sigma2 :
						res.transpose().multiply(res).getValueAt(0, 0) / upsilon;
			return new Model(this, betaHat, invXtX, sigma2Hat, upsilon);
		}
		
	}
	
	static class Model {
		final Matrix betaHat;
		final SymmetricMatrix invXtX;
		final double sigma2Hat;
		final double sigmaHat;
		Matrix invXtXChol;
		final int upsilon;
		final Sample sample;
		final Matrix residuals;
		final Map<Integer, List<Integer>> sampleIndexMap;
		
		Model(Sample sample, Matrix betaHat, SymmetricMatrix invXtX, double sigma2Hat, int upsilon) {
			this.sample = sample;
			this.betaHat = betaHat;
			this.invXtX = invXtX; 
			this.sigma2Hat = sigma2Hat;
			this.upsilon = upsilon;
			this.sigmaHat = Math.sqrt(sigma2Hat);
			residuals = sample.getVectorY().subtract(sample.getMatrixX().multiply(betaHat));
			sampleIndexMap = new HashMap<Integer, List<Integer>>();
		}
		
		private Matrix getOmegaChol() {
			if (invXtXChol == null) {
				invXtXChol = invXtX.getLowerCholTriangle();
			}
			return invXtXChol;
		}
		
		private Matrix sampleResiduals(int n) {
			if (!sampleIndexMap.containsKey(n)) {
				List<Integer> sampleIndex = new ArrayList<Integer>();
				for (int i = 0; i < n; i++)
					sampleIndex.add(i);
				sampleIndexMap.put(n, sampleIndex);
			}
			
			List<Integer> index = SamplingUtility.getSample(sampleIndexMap.get(n), n, true);
			Matrix sampledResiduals = residuals.getSubMatrix(index, null, false); // we do not sort the indices
			return sampledResiduals;
		}
		
		Matrix getRandomDeviate(List<Double> xValues, Transformation t, boolean useNonParametricMethod) {
			Matrix xMat = createMatrixX(xValues);
			Matrix betaDeviates = getOmegaChol().multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN)).scalarMultiply(sigmaHat);
			Matrix betaHat_b = betaHat.add(betaDeviates);
			
			Matrix epsilon_b = useNonParametricMethod ? 
					sampleResiduals(xMat.m_iRows) :
						StatisticalUtility.drawRandomVector(xMat.m_iRows, Type.GAUSSIAN).scalarMultiply(sigmaHat);
			
			Matrix result = xMat.multiply(betaHat_b).add(epsilon_b);
			return convertPredictionToOriginalScale(result, t);
		}
		
		private Matrix convertPredictionToOriginalScale(Matrix transScalePred, Transformation t) {
			switch(t) {
			case Sqr:
				return transScalePred.elementWisePower(2d);
			case Exp:
				return transScalePred.expMatrix();
			case Log:
				return transScalePred.logMatrix();
			case Sqrt:
				return transScalePred.elementWisePower(.5);
			default:
				throw new InvalidParameterException("This transformation " + t.name() + " is not supported!");
			}
		}

		Matrix createMatrixX(List<Double> xValues) {
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

		private Matrix getBeauchampAndOlsonEstimator(List<Double> xValues) {
			Matrix xMat = createMatrixX(xValues);
			Matrix xBeta = getXBeta(xValues);
			
			Matrix output = new Matrix(xValues.size(), 1);
			
			for (int ii = 0; ii < xValues.size(); ii++) {
				double term1 = Math.exp(xBeta.getValueAt(ii, 0) + sigma2Hat * .5);
				Matrix xMat_ii = xMat.getSubMatrix(ii, ii, 0, xMat.m_iCols - 1);
				int n = upsilon + xMat_ii.m_iCols;
				Matrix xInvXtXxT = xMat_ii.multiply(invXtX).multiply(xMat_ii.transpose());
				double phi = xInvXtXxT.getValueAt(0, 0) * n;
				double sigma2hat2 = sigma2Hat * sigma2Hat;
				double term2 = 1 - sigma2Hat * (2*phi + sigma2Hat)/(4*n) + 
						sigma2hat2 * (sigma2hat2 + 2 * (16d/3 + 2 * phi) * sigma2Hat + 4 * phi * phi + 16 * phi) / (32 * n * n);
				output.setValueAt(ii, 0, term1 * term2);
			}
			return output;
		}

		private Matrix getBeauchampAndOlsonEstimatorVariance(List<Double> xValues) {
			Matrix xMat = createMatrixX(xValues);
			Matrix xBeta = getXBeta(xValues);
			
			Matrix output = new Matrix(xValues.size(), 1);
			
			for (int ii = 0; ii < xValues.size(); ii++) {
				double term1 = Math.exp(2*xBeta.getValueAt(ii, 0) + sigma2Hat);
				Matrix xMat_ii = xMat.getSubMatrix(ii, ii, 0, xMat.m_iCols - 1);
				int n = upsilon + xMat_ii.m_iCols;
				Matrix xInvXtXxT = xMat_ii.multiply(invXtX).multiply(xMat_ii.transpose());
				double phi = xInvXtXxT.getValueAt(0, 0) * n;
				double sigma2hat2 = sigma2Hat * sigma2Hat;
				double term2 = Math.exp(sigma2Hat) * (1 - 2 * sigma2Hat * (phi + 2*sigma2Hat)/n + 
						2* sigma2hat2 * (4*sigma2hat2 + 2 * (16d/3 + 2 * phi) * sigma2Hat + phi * phi + 4 * phi) / (n * n));
				double term3 = 1 - sigma2Hat/n * (sigma2Hat + 2 * phi) + sigma2hat2 / (2*n*n) * (sigma2hat2 + (16d/3 + 4*phi) * sigma2Hat + 4 * phi * phi + 8 * phi);
				output.setValueAt(ii, 0, term1 * (term2 - term3));
			}
			return output;
		}

		
		private Matrix getBaskervilleEstimator(List<Double> xValues) {
			Matrix xBeta = getXBeta(xValues);
			return xBeta.scalarAdd(sigma2Hat * .5).expMatrix();
		}
		
		private Matrix getBaskervilleEstimatorVariance(List<Double> xValues) {
			Matrix xBeta = getXBeta(xValues).scalarMultiply(2d).scalarAdd(sigma2Hat).expMatrix();
			return xBeta.scalarMultiply(Math.exp(sigma2Hat) - 1);
		}

//		private Matrix getSmearingEstimate(List<Double> xValues, Transformation t) {
//			Matrix xBeta = getXBeta(xValues);
//			Matrix res_hat = sample.getVectorY().subtract(sample.getMatrixX().multiply(betaHat));
//			Matrix result = null;
//			for (int i = 0; i < res_hat.m_iRows; i++) {
//				result = result == null ?
//						convertPredictionToOriginalScale(xBeta.scalarAdd(res_hat.getValueAt(i, 0)), t) :
//							result.add(convertPredictionToOriginalScale(xBeta.scalarAdd(res_hat.getValueAt(i, 0)), t));	
//			}
//			return result.scalarMultiply(1d / res_hat.m_iRows);
//		}
		
		ComplexMatrix getComplexRandomDeviate(List<Double> xValues, Transformation t, boolean useNonParametricMethod) {
			double chiSquareDeviate = StatisticalUtility.getRandom().nextChiSquare(upsilon);
			ComplexNumber sigma2Hat_b = new ComplexNumber(sigma2Hat, sigma2Hat * (chiSquareDeviate - upsilon) / upsilon);
			ComplexNumber sigmaHat_b = sigma2Hat_b.sqrt();

			Matrix xMat = createMatrixX(xValues);

			Matrix betaDeviatesMat = getOmegaChol().multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN));
			
			ComplexNumber[] cnArray = new ComplexNumber[xValues.size()];
			for (int ii = 0; ii < xValues.size(); ii++) {
				Matrix xMat_ii = xMat.getSubMatrix(ii, ii, 0, xMat.m_iCols - 1);
				ComplexNumber betaDeviates = new ComplexNumber(0, xMat_ii.multiply(betaDeviatesMat).getValueAt(0, 0));
				betaDeviates = betaDeviates.multiply(sigmaHat_b);

				Number epsilon_b = useNonParametricMethod ?
						sampleResiduals(1).getValueAt(0, 0) :
							sigmaHat_b.multiply(StatisticalUtility.getRandom().nextGaussian());
				double xBetaHat = xMat_ii.multiply(betaHat).getValueAt(0, 0); 
				ComplexNumber mean = betaDeviates.add(xBetaHat);
				ComplexNumber result = mean.add(epsilon_b);
				switch(t) {
				case Sqr:
					cnArray[ii] = result.square();
					break;
				case Exp:
					cnArray[ii] = result.exp();
					break;
				case Log:
					cnArray[ii] = result.log();
					break;
				case Sqrt:
					cnArray[ii] = result.sqrt();
					break;
				default:
					throw new InvalidParameterException("The transformation " + t.name() + " is not supported!");
				}
			}
			return new ComplexMatrix(cnArray);
		}

		private Matrix getGregoireEstimator(List<Double> xValues) {
			Matrix xMat = createMatrixX(xValues);
			Matrix varXBeta = xMat.multiply(invXtX).multiply(xMat.transpose()).scalarMultiply(sigma2Hat).diagonalVector();
			Matrix pred = getXBeta(xValues);
			
			return pred.elementWisePower(2d).subtract(varXBeta).scalarAdd(sigma2Hat);
		}
	}
	
	static enum Transformation {Sqr, Exp, Log, Sqrt}
	
	

	private static void doRun(Transformation t, int sampleSize, int nbRealizations) throws IOException {
		System.out.println("Simulating " + t.name() + " with sample size n = " + sampleSize);
		int nbInnerReal = 10000;
		String filename = ObjectUtility.getPackagePath(ComplexNumberSimpleCaseStudy.class) + "caseStudy_" + t.name() + sampleSize +".csv";
		filename = filename.replace("bin/", "");
		CSVWriter writer = new CSVWriter(new File(filename), false);
		List<FormatField> fields = new ArrayList<FormatField>();
		fields.add(new CSVField("RealID"));
		fields.add(new CSVField("b0"));
		fields.add(new CSVField("b1"));
		fields.add(new CSVField("sigma2"));
		fields.add(new CSVField("x"));
		fields.add(new CSVField("MC_OLD"));
		fields.add(new CSVField("MC_OLDVar"));
		fields.add(new CSVField("MC_NEW"));
		fields.add(new CSVField("MC_NEW_imag"));
		fields.add(new CSVField("MC_NEWVar_real"));
		fields.add(new CSVField("MC_NEWVar_imag"));
		if (t == Transformation.Exp) {
			fields.add(new CSVField("Beauchamp"));
			fields.add(new CSVField("BeauchampVar"));
			fields.add(new CSVField("Baskerville"));
			fields.add(new CSVField("BaskervilleVar"));
		} else if (t == Transformation.Sqr) {
			fields.add(new CSVField("Gregoire"));
		}
		writer.setFields(fields);
		
		Matrix trueBeta = new Matrix(2,1);
		double trueVariance;
		if (t == Transformation.Exp) {
			trueBeta.setValueAt(0, 0, 2);
			trueBeta.setValueAt(1, 0, 0.25);
			trueVariance = 1d;
		} else if (t == Transformation.Sqr) {  
			trueBeta.setValueAt(0, 0, 2);
			trueBeta.setValueAt(1, 0, 0.25);
			trueVariance = 2d;
		} else {
			trueBeta.setValueAt(0, 0, 2.5);
			trueBeta.setValueAt(1, 0, 4);
			trueVariance = 12d;
		}
		List<Double> xValues = new ArrayList<Double>();
		for (double i = 3; i <= 10; i++) {
			xValues.add(i);
		}
		for (int real = 1; real <= nbRealizations; real++) {
			if (real % 1000 == 0) {
				System.out.println("Running realization " + real);
			}
			Sample s = Sample.createSample(trueBeta, trueVariance, sampleSize);
			Model m = s.getModel(null); // null means we are relying on the estimated variance
			Matrix beauchampAndOlsonEstimator = null;
			Matrix beauchampAndOlsonEstimatorVar = null;
			Matrix baskervilleEstimator = null;
			Matrix baskervilleEstimatorVar = null;
			Matrix gregoireEstimator = null;
			if (t == Transformation.Exp) {
				beauchampAndOlsonEstimator = m.getBeauchampAndOlsonEstimator(xValues);
				beauchampAndOlsonEstimatorVar = m.getBeauchampAndOlsonEstimatorVariance(xValues);
				baskervilleEstimator = m.getBaskervilleEstimator(xValues);
				baskervilleEstimatorVar = m.getBaskervilleEstimatorVariance(xValues);
			} else if (t == Transformation.Sqr) {
				gregoireEstimator = m.getGregoireEstimator(xValues);
			}
						
			MonteCarloEstimate mcEstimator = new MonteCarloEstimate();
			ComplexMonteCarloEstimate cmcEstimator = new ComplexMonteCarloEstimate();

			for (int innerReal = 0; innerReal < nbInnerReal; innerReal++) {
				mcEstimator.addRealization(m.getRandomDeviate(xValues, t, false)); // false: regular parametric method
				cmcEstimator.addRealization(m.getComplexRandomDeviate(xValues, t, false)); // false : regular parametric method
			}
			
			Matrix mcMean = mcEstimator.getMean();
			Matrix mcVar = mcEstimator.getVariance();
			
			ComplexMatrix cmcMean = cmcEstimator.getMean();
			Matrix cmcRealPart = cmcEstimator.getVarianceRealPart();
			Matrix cmcImagPart = cmcEstimator.getVarianceImaginaryPart();
			
			for (int ii = 0; ii < xValues.size(); ii++) {
				Object[] record = new Object[fields.size()];
				record[0] = real;
				record[1] = m.betaHat.getValueAt(0, 0);
				record[2] = m.betaHat.getValueAt(1, 0);
				record[3] = m.sigma2Hat;
				record[4] = xValues.get(ii);
				record[5] = mcMean.getValueAt(ii, 0);
				record[6] = mcVar.getValueAt(ii, ii);
				record[7] = cmcMean.getValueAt(ii, 0).realPart;
				record[8] = cmcMean.getValueAt(ii, 0).imaginaryPart;
				record[9] = cmcRealPart.getValueAt(ii, 0);
				record[10] = cmcImagPart.getValueAt(ii, 0);
				if (t == Transformation.Exp) {
					record[11] = beauchampAndOlsonEstimator.getValueAt(ii, 0);
					record[12] = beauchampAndOlsonEstimatorVar.getValueAt(ii, 0);
					record[13] = baskervilleEstimator.getValueAt(ii, 0);
					record[14] = baskervilleEstimatorVar.getValueAt(ii, 0);
				} else  if (t == Transformation.Sqr) {
					record[9] = gregoireEstimator.getValueAt(ii, 0);
				}
				writer.addRecord(record);
			}
		}
		writer.close();

	}
	
//	private static void runSimulationChiSquare(int nReal) {
//		int upsilon = 50;
//		double A_j = 0.8;
//		double sigma2_hat = 1;
//		double realPart = 0;
//		double imagPart = 0;
//		for (int i = 0; i < nReal; i++) {
//			ComplexNumber exponent = new ComplexNumber(0, A_j * sigma2_hat / upsilon * (StatisticalUtility.getRandom().nextChiSquare(upsilon) - upsilon));
//			ComplexNumber result = exponent.exp();
//			realPart += result.realPart;
//			imagPart += result.imaginaryPart;
//		}
//		System.out.println("Real part = " + (realPart/nReal));
//		System.out.println("Imaginary part = " + (imagPart/nReal));
//		double approx = Math.exp(- A_j * A_j * sigma2_hat * sigma2_hat / upsilon); 
//		System.out.println("Approximation = " + approx);
//	}

	
	private static void runSimulationChiSquareToSeries(int nReal) {
		int upsilon = 100;
		double sigma2 = 2d;
		double mean = 0d;
		for (int i = 0; i < nReal; i++) {
			double sigma2_hat = sigma2 * StatisticalUtility.getRandom().nextChiSquare(upsilon) / upsilon;
			mean += Math.exp(sigma2_hat/2);
		}
		System.out.println("Monte Carlo = " + (mean/nReal));

		double series = 0d;
		for (int k = 0; k < 20; k++) {
			series += Math.pow(sigma2/2, k) / MathUtility.Factorial(k) * Math.pow(2,k) * GammaUtility.gamma(k + 0.5 * upsilon) / (Math.pow(upsilon, k) * GammaUtility.gamma(0.5 * upsilon));
			System.out.println("Series (k = " + k + ") = " + series);
		}
		System.out.println("With sigma2 = " + (Math.exp(sigma2/2)));
	}

	public static void main(String[] arg) throws IOException {
//		runSimulationChiSquareToSeries(1000000);
//		runSimulationChiSquare(1000000);
		doRun(Transformation.Exp, 50, 10000);
//		doRun(Transformation.Exp, 100);
//		doRun(Transformation.Exp, 200);
//		doRun(Transformation.Sqr, 50);
//		doRun(Transformation.Sqr, 100);
//		doRun(Transformation.Sqr, 200);
	}	
}
