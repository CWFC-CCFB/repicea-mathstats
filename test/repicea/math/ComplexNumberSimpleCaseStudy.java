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
import java.util.List;

import repicea.io.FormatField;
import repicea.io.javacsv.CSVField;
import repicea.io.javacsv.CSVWriter;
import repicea.stats.Distribution;
import repicea.stats.StatisticalUtility;
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
		
		double getRandomDeviate(double x) {
			Matrix xMat = new Matrix(1,2);
			xMat.setValueAt(0, 0, 1d);
			xMat.setValueAt(0, 1, x);
			Matrix betaDeviates = getOmegaChol().multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN)).scalarMultiply(sigmaHat);
			
			Matrix betaHat_b = betaHat.add(betaDeviates);
			double epsilon_b = StatisticalUtility.getRandom().nextGaussian() * sigmaHat;
			return betaHat_b.getValueAt(0, 0) + betaHat_b.getValueAt(1, 0) * x + epsilon_b;
		}
		
		ComplexNumber getComplexRandomDeviate(double x) {
			Matrix xMat = new Matrix(1,2);
			xMat.setValueAt(0, 0, 1d);
			xMat.setValueAt(0, 1, x);
			
			double chiSquareDeviate = StatisticalUtility.getRandom().nextChiSquare(upsilon);
			ComplexNumber sigma2Hat_b = new ComplexNumber(sigma2Hat, sigma2Hat * (chiSquareDeviate - upsilon) / upsilon);
			ComplexNumber sigmaHat_b = sigma2Hat_b.sqrt();
			
			Matrix betaDeviatesMat = getOmegaChol().multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN));
			ComplexNumber betaDeviates = new ComplexNumber(0, xMat.multiply(betaDeviatesMat).getValueAt(0, 0));
			betaDeviates = betaDeviates.multiply(sigmaHat_b);
			
			ComplexNumber epsilon_b = sigmaHat_b.multiply(StatisticalUtility.getRandom().nextGaussian());
			
			ComplexNumber mean = betaDeviates.add(betaHat.getValueAt(0, 0) + betaHat.getValueAt(1, 0) * x);
			return mean.add(epsilon_b);
		}
	}
	
	static enum Transformation {Sqr, Exp, Log, Sqrt}
	
	private static Number computeValue(Transformation t, Number originalValue) {
		switch(t) {
		case Sqr:
			return originalValue instanceof ComplexNumber ?
					((ComplexNumber) originalValue).square() :
						originalValue.doubleValue() * originalValue.doubleValue();
		case Exp:
			return originalValue instanceof ComplexNumber ?
					((ComplexNumber) originalValue).exp() :
						Math.exp(originalValue.doubleValue());
		case Log:
			return originalValue instanceof ComplexNumber ?
					((ComplexNumber) originalValue).log() :
						Math.log(originalValue.doubleValue());
		case Sqrt:
			return originalValue instanceof ComplexNumber ?
					((ComplexNumber) originalValue).sqrt() :
						Math.sqrt(originalValue.doubleValue());
		default:
			throw new InvalidParameterException();
		}
	}
	
	
	private static void doRun(Transformation t, int sampleSize) throws IOException {
		System.out.println("Simulating " + t.name() + " with sample size n = " + sampleSize);
		int nbRealizations = 10000;
		int nbInnerReal = 10000;
		String filename = ObjectUtility.getPackagePath(ComplexNumberSimpleCaseStudy.class) + "caseStudy_" + t.name() + sampleSize +".csv";
		filename = filename.replace("bin/", "");
		CSVWriter writer = new CSVWriter(new File(filename), false);
		List<FormatField> fields = new ArrayList<FormatField>();
		fields.add(new CSVField("RealID"));
		fields.add(new CSVField("b0"));
		fields.add(new CSVField("b1"));
		fields.add(new CSVField("sigma2"));
		fields.add(new CSVField("yOLD"));
		fields.add(new CSVField("yNEW"));
		fields.add(new CSVField("yNEW_imag"));
		writer.setFields(fields);
		
		Matrix trueBeta = new Matrix(2,1);
		double trueVariance;
		if (t == Transformation.Sqr || t == Transformation.Exp) {
			trueBeta.setValueAt(0, 0, 2);
			trueBeta.setValueAt(1, 0, 0.25);
			trueVariance = 4d;
		} else {
			trueBeta.setValueAt(0, 0, 2.5);
			trueBeta.setValueAt(1, 0, 4);
			trueVariance = 12d;
		}
		
		for (int real = 1; real <= nbRealizations; real++) {
			if (real % 1000 == 0) {
				System.out.println("Running realization " + real);
			}
			Sample s = Sample.createSample(trueBeta, trueVariance, sampleSize);
			Model m = s.getModel(null); // null means we are relying on the estimated variance
			double sumOLD = 0;
			ComplexNumber sumNEW = new ComplexNumber(0, 0);
			for (int innerReal = 0; innerReal < nbInnerReal; innerReal++) {
				double n1 = m.getRandomDeviate(5);
				sumOLD += computeValue(t, n1).doubleValue();
				ComplexNumber n2 = m.getComplexRandomDeviate(5);
				sumNEW = sumNEW.add(computeValue(t, n2));
			}
			sumOLD /= nbInnerReal;
			sumNEW = sumNEW.multiply(1d / nbInnerReal);
			Object[] record = new Object[7];
			record[0] = real;
			record[1] = m.betaHat.getValueAt(0, 0);
			record[2] = m.betaHat.getValueAt(1, 0);
			record[3] = m.sigma2Hat;
			record[4] = sumOLD;
			record[5] = sumNEW.doubleValue();
			record[6] = sumNEW.imaginaryPart;
			writer.addRecord(record);
		}
		writer.close();

	}
	
	public static void main(String[] arg) throws IOException {
//		doRun(Transformation.Sqr, 50);
//		doRun(Transformation.Sqr, 100);
//		doRun(Transformation.Sqr, 200);
//		doRun(Transformation.Exp, 50);
//		doRun(Transformation.Exp, 100);
//		doRun(Transformation.Exp, 200);
		doRun(Transformation.Log, 50);
		doRun(Transformation.Log, 100);
		doRun(Transformation.Log, 200);
		doRun(Transformation.Sqrt, 50);
		doRun(Transformation.Sqrt, 100);
		doRun(Transformation.Sqrt, 200);
	}	
}
