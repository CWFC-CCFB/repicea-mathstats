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
import java.util.ArrayList;
import java.util.List;

import repicea.io.FormatField;
import repicea.io.javacsv.CSVField;
import repicea.io.javacsv.CSVWriter;
import repicea.stats.Distribution;
import repicea.stats.StatisticalUtility;
import repicea.util.ObjectUtility;

/**
 * The ComplexNumberDiffEquationCaseStudy class implements a simulation study on
 * differential growth equation. 
 * @author Mathieu Fortin - Oct 2023
 */
public class ComplexNumberDiffEquationCaseStudy {

	
	static class PopulationUnit {
		final double x;
		final double y;
		final int id;
		PopulationUnit(Matrix trueBeta, double trueStd, int id) {
			this.id = id;
			x = StatisticalUtility.getRandom().nextDouble() * 140;
			y = trueBeta.getValueAt(0, 0) + x * trueBeta.getValueAt(1, 0) + x * x * trueBeta.getValueAt(2, 0) + StatisticalUtility.getRandom().nextGaussian() * trueStd;
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
			Matrix xMat = new Matrix(size(), 3);
			for (int i = 0; i < size(); i++) {
				xMat.setValueAt(i, 0, 1d);
				xMat.setValueAt(i, 1, get(i).x);
				xMat.setValueAt(i, 2, get(i).x * get(i).x);
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
			return new Model(betaHat, invXtX, sigma2Hat, upsilon);
		}
		
	}
	
	static class Model {
		final Matrix betaHat;
		final SymmetricMatrix invXtX;
		final double sigma2Hat;
		final double sigmaHat;
		Matrix invXtXChol;
		final int upsilon;
		
		Model(Matrix betaHat, SymmetricMatrix invXtX, double sigma2Hat, int upsilon) {
			this.betaHat = betaHat;
			this.invXtX = invXtX; 
			this.sigma2Hat = sigma2Hat;
			this.upsilon = upsilon;
			this.sigmaHat = Math.sqrt(sigma2Hat);
		}
		
		private Matrix getOmegaChol() {
			if (invXtXChol == null) {
				invXtXChol = invXtX.getLowerCholTriangle();
			}
			return invXtXChol;
		}
		
		Matrix getRandomDeviate(double x, int nbCycles) {
			Matrix betaDeviates = getOmegaChol().multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN)).scalarMultiply(sigmaHat);
			
			Matrix betaHat_b = betaHat.add(betaDeviates);
			Matrix output = new Matrix(nbCycles, 1);
			for (int i = 0; i < nbCycles; i++) {
				double epsilon_b = StatisticalUtility.getRandom().nextGaussian() * sigmaHat;
				x = betaHat_b.getValueAt(0, 0) + betaHat_b.getValueAt(1, 0) * x + betaHat_b.getValueAt(2, 0) * x * x + epsilon_b;
				output.setValueAt(i, 0, x);
			}
			return output;
		}

		static Matrix getTrueRealization(double x, int nbCycles, Matrix beta, double sigma2) {
			double sigma = Math.sqrt(sigma2);
			Matrix output = new Matrix(nbCycles, 1);
			for (int i = 0; i < nbCycles; i++) {
				double epsilon = StatisticalUtility.getRandom().nextGaussian() * sigma;
				x = beta.getValueAt(0, 0) + beta.getValueAt(1, 0) * x + beta.getValueAt(2, 0) * x * x + epsilon;
				output.setValueAt(i, 0, x);
			}
			return output;
		}

		ComplexNumber[] getComplexRandomDeviate(double x, int nbCycles) {
			Matrix xMat = new Matrix(1,3);
			xMat.setValueAt(0, 0, 1d);
			xMat.setValueAt(0, 1, x);
			xMat.setValueAt(0, 2, x*x);
			
			double chiSquareDeviate = StatisticalUtility.getRandom().nextChiSquare(upsilon);
			ComplexNumber sigma2Hat_b = new ComplexNumber(sigma2Hat, sigma2Hat * (chiSquareDeviate - upsilon) / upsilon);
			ComplexNumber sigmaHat_b = sigma2Hat_b.sqrt();
			
			Matrix betaDeviatesMat = getOmegaChol().multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN));
			ComplexNumber betaDeviates = new ComplexNumber(0, xMat.multiply(betaDeviatesMat).getValueAt(0, 0));
			betaDeviates = betaDeviates.multiply(sigmaHat_b);
			Number xNum = x;
			ComplexNumber[] output = new ComplexNumber[nbCycles];
			for (int i = 0; i < nbCycles; i++) {
				ComplexNumber epsilon_b = sigmaHat_b.multiply(StatisticalUtility.getRandom().nextGaussian());
				Number xBeta1 = xNum instanceof ComplexNumber ?
						((ComplexNumber) xNum).multiply(betaHat.getValueAt(1, 0)) :
							xNum.doubleValue() * betaHat.getValueAt(1, 0);
				Number xBeta2 = xNum instanceof ComplexNumber ?
						((ComplexNumber) xNum).square().multiply(betaHat.getValueAt(2, 0)) :
							xNum.doubleValue() * xNum.doubleValue() * betaHat.getValueAt(2, 0);
				ComplexNumber mean = betaDeviates.add(betaHat.getValueAt(0, 0));
				
				mean = mean.add(xBeta1);
				mean = mean.add(xBeta2);
				mean = mean.add(epsilon_b);
				output[i] = mean;
				xNum = mean;
				
			}
			return output;
		}
	}
	
	public static void main(String[] arg) throws IOException {
		int nbRealizations = 1000;
		int nbInnerReal = 10000;
		int sampleSize = 50;
		int nbCycles = 10;
		String filename = ObjectUtility.getPackagePath(ComplexNumberDiffEquationCaseStudy.class) + "caseStudyDiffEq_" + sampleSize +".csv";
		filename = filename.replace("bin/", "");
		CSVWriter writer = new CSVWriter(new File(filename), false);
		List<FormatField> fields = new ArrayList<FormatField>();
		fields.add(new CSVField("RealID"));
		fields.add(new CSVField("b0"));
		fields.add(new CSVField("b1"));
		fields.add(new CSVField("b2"));
		fields.add(new CSVField("sigma2"));
		for (int j = 1; j <= nbCycles; j++) {
			fields.add(new CSVField("yOLD" + j));
			fields.add(new CSVField("yNEW" + j));
			fields.add(new CSVField("yNEW_imag" + j));
		}
		writer.setFields(fields);
		
		Matrix trueBeta = new Matrix(3,1);
		trueBeta.setValueAt(0, 0, 15);
		trueBeta.setValueAt(1, 0, 1.4);
		trueBeta.setValueAt(2, 0, -0.004);
		double trueVariance = 16d;
		
		for (int real = 0; real < nbRealizations; real++) {
			System.out.println("Running realization " + (real + 1));
			Sample s = Sample.createSample(trueBeta, trueVariance, sampleSize);
			Model m = s.getModel(null); // null means we are relying on the estimated variance
			Matrix sumOLD = null;
			ComplexNumber[] sumNEW = null;
			for (int innerReal = 0; innerReal < nbInnerReal; innerReal++) {
				sumOLD = sumOLD == null ?
						m.getRandomDeviate(0, nbCycles) :
							sumOLD.add(m.getRandomDeviate(0, nbCycles));
				ComplexNumber[] n2 = m.getComplexRandomDeviate(0, nbCycles);
				if (sumNEW == null) {
					sumNEW = n2;
				} else {
					for (int j = 0; j < sumNEW.length; j++) {
						sumNEW[j] = sumNEW[j].add(n2[j]);
					}
				}
						
			}
			sumOLD = sumOLD.scalarMultiply(1d / nbInnerReal);
			for (int j = 0; j < sumNEW.length; j++) {
				sumNEW[j] = sumNEW[j].multiply(1d / nbInnerReal);
			}
			Object[] record = new Object[nbCycles * 3 + 5];
			record[0] = real;
			record[1] = m.betaHat.getValueAt(0, 0);
			record[2] = m.betaHat.getValueAt(1, 0);
			record[3] = m.betaHat.getValueAt(2, 0);
			record[4] = m.sigma2Hat;
			for (int j = 0; j < nbCycles; j++) {
				record[j*3 + 5] = sumOLD.getValueAt(j, 0);
				record[j*3 + 6] = sumNEW[j].doubleValue(); 
				record[j*3 + 7] = sumNEW[j].imaginaryPart; 
			}
			writer.addRecord(record);
		}
		Matrix sum = null;
		for (int innerReal = 0; innerReal < 100000; innerReal++) {
			Matrix real = Model.getTrueRealization(0, nbCycles, trueBeta, trueVariance);
			sum = sum == null ?
					real :
						sum.add(real);
		}
		sum = sum.scalarMultiply(1d / 100000);
		Object[] record = new Object[nbCycles * 3 + 5];
		record[0] = 100000;
		record[1] = trueBeta.getValueAt(0, 0);
		record[2] = trueBeta.getValueAt(1, 0);
		record[3] = trueBeta.getValueAt(2, 0);
		record[4] = trueVariance;
		for (int j = 0; j < nbCycles; j++) {
			record[j*3 + 5] = sum.getValueAt(j, 0);
			record[j*3 + 6] = 0; 
			record[j*3 + 7] = 0; 
		}
		writer.addRecord(record);
		
		writer.close();
	}	
}
