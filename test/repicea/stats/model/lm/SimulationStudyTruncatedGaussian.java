/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2023 His Majesty the King in Right of Canada
 * Author: Mathieu Fortin - Canadian Wood Fibre Centre, CFS
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
package repicea.stats.model.lm;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import repicea.io.FormatField;
import repicea.io.javacsv.CSVField;
import repicea.io.javacsv.CSVWriter;
import repicea.math.Matrix;
import repicea.stats.REpiceaRandom;
import repicea.stats.StatisticalUtility;
import repicea.stats.data.DataSet;
import repicea.util.ObjectUtility;

class SimulationStudyTruncatedGaussian {

	static final REpiceaRandom RANDOM = StatisticalUtility.getRandom();
	static final double b0 = -1d;
	static final double b1 = .5;
	
	class ObservationalUnit {
		private final double x;
		private final double w;
		private final double y;
		
		ObservationalUnit() {
			x = RANDOM.nextDouble() * (deltaX) + xLower;
			double xBeta = b0 + b1 * x;
			double potentialWValue;
			while ((potentialWValue = xBeta + RANDOM.nextGaussian() * sigma) < 0) {};
			w = potentialWValue;
			y = Math.exp(w) - 1;
		}
		
		ObservationalUnit(double x) {
			this.x = x;
			double xBeta = b0 + b1 * x;
			double potentialWValue;
			while ((potentialWValue = xBeta + RANDOM.nextGaussian() * sigma) < 0) {};
			w = potentialWValue;
			y = Math.exp(w) - 1;
		}
		
		
		
		Object[] convertToObservation() {
			Object[] record = new Object[2];
			record[0] = w;
			record[1] = x;
			return record;
		}
	}


	static final double sigma2 = 2d;
	static final double sigma = Math.sqrt(sigma2);
	
	private List<ObservationalUnit> samplePopulation(int size) {
		List<ObservationalUnit> units = new ArrayList<ObservationalUnit>();
		while(units.size() < size) {
			units.add(new ObservationalUnit());
		}
		return units;
	}
	
	private DataSet createDataSet(int size) {
		DataSet ds = new DataSet(Arrays.asList(new String[] {"w", "x"}));
		for (ObservationalUnit u : samplePopulation(size)) {
			ds.addObservation(u.convertToObservation());
		}
		ds.indexFieldType();
		return ds;
	}
	
	final int nbRealizations;
	final int sampleSize;
	final double xLower;
	final double xUpper;
	final String filename;
	final double deltaX;
	
	SimulationStudyTruncatedGaussian(int nbRealizations, int sampleSize, double xLower, double xUpper, String filename) {
		this.nbRealizations = nbRealizations;
		this.sampleSize = sampleSize;
		this.xLower = xLower;
		this.xUpper = xUpper;
		this.deltaX = xUpper - xLower;
		this.filename = filename;
	}
	
	public void run() throws IOException {
		CSVWriter writer = new CSVWriter(new File(filename), false);
		List<FormatField> fields = new ArrayList<FormatField>();
		fields.add(new CSVField("b0_trad"));
		fields.add(new CSVField("b1_trad"));
		fields.add(new CSVField("sigma2_trad"));
		fields.add(new CSVField("b0_trad_var"));
		fields.add(new CSVField("b1_trad_var"));
		fields.add(new CSVField("b0_new"));
		fields.add(new CSVField("b1_new"));
		fields.add(new CSVField("sigma2_new"));
		fields.add(new CSVField("b0_new_var"));
		fields.add(new CSVField("b1_new_var"));
		fields.add(new CSVField("sigma2_new_var"));
		writer.setFields(fields);
		
		for (int i = 0; i < nbRealizations; i++) {
			if ((i+1) % 100 == 0) {
				System.out.println("Processing realization " + (i+1));
			}
			DataSet ds = createDataSet(sampleSize);
			LinearModel lm = new LinearModel(ds, "w ~ x");
			lm.doEstimation();
			Matrix p = lm.getParameters();
			Matrix vTrad = lm.getEstimator().getParameterEstimates().getVariance().diagonalVector();
			Matrix parmsTrad = new Matrix(p.m_iRows + 1, 1);
			parmsTrad.setSubMatrix(p, 0, 0);
			parmsTrad.setValueAt(parmsTrad.m_iRows - 1, 0, lm.getResidualVariance());
			LinearModelWithTruncatedGaussianErrorTerm truncatedModel = new LinearModelWithTruncatedGaussianErrorTerm(ds, "w ~ x", parmsTrad, 0);
			truncatedModel.doEstimation();
			if (truncatedModel.getEstimator().isConvergenceAchieved()) {
				Matrix truncatedParms = truncatedModel.getParameters();
				Matrix truncatedVar = truncatedModel.getEstimator().getParameterEstimates().getVariance().diagonalVector();
				Matrix output = parmsTrad.matrixStack(vTrad, true).matrixStack(truncatedParms, true).matrixStack(truncatedVar, true);
				Object[] record = new Object[output.m_iRows];
				for (int j = 0; j < output.m_iRows; j++) 
					record[j] = output.getValueAt(j, 0);
				writer.addRecord(record);
			}
		}
		writer.close();
	}
	
	private void produceMeanValues(String filename) throws IOException {
		CSVWriter writer = new CSVWriter(new File(filename), false);
		List<FormatField> fields = new ArrayList<FormatField>();
		fields.add(new CSVField("x"));
		fields.add(new CSVField("w"));
		fields.add(new CSVField("y"));
		writer.setFields(fields);
		
		for (double x = xLower; x <= xUpper; x+=.5) {
			System.out.println("Processing x = " + x);
			List<ObservationalUnit> units = new ArrayList<ObservationalUnit>();
			while(units.size() < 1000000) {
				units.add(new ObservationalUnit(x));
			}
			double meanW = 0;
			double meanY = 0;
			for (ObservationalUnit u : units) {
				meanY += u.y;
				meanW += u.w;
			}
			meanY /= units.size();
			meanW /= units.size();
			Object[] record = new Object[3];
			record[0] = x;
			record[1] = meanW;
			record[2] = meanY;
			writer.addRecord(record);
		}
		writer.close();
	}
		
	
	

	public static void main(String[] args) throws IOException {
		String filename = ObjectUtility.getPackagePath(SimulationStudyTruncatedGaussian.class) + "meanValues.csv";
		filename = filename.replace("/bin", "");
		System.out.println("Saving mean values in " + filename);
		SimulationStudyTruncatedGaussian sim = new SimulationStudyTruncatedGaussian(10000, 500, 3, 10, filename);
		sim.produceMeanValues(filename);

//		filename = ObjectUtility.getPackagePath(SimulationStudyTruncatedGaussian.class) + "Simulation100.csv";
//		filename = filename.replace("/bin", "");
//		System.out.println("Saving simulation result in " + filename);
//		sim = new SimulationStudyTruncatedGaussian(10000, 100, 3, 10, filename);
//		sim.run();
//
//		filename = ObjectUtility.getPackagePath(SimulationStudyTruncatedGaussian.class) + "Simulation500.csv";
//		filename = filename.replace("/bin", "");
//		System.out.println("Saving simulation result in " + filename);
//		sim = new SimulationStudyTruncatedGaussian(10000, 500, 3, 10, filename);
//		sim.run();
	}
	
	
}
