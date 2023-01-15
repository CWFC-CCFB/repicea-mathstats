/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2022 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.model.glm.measerr.simulation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;

import repicea.io.FormatField;
import repicea.io.javacsv.CSVField;
import repicea.io.javacsv.CSVWriter;
import repicea.math.Matrix;
import repicea.serial.xml.XmlDeserializer;
import repicea.serial.xml.XmlSerializer;
import repicea.stats.estimators.MaximumLikelihoodEstimator;
import repicea.stats.model.AbstractStatisticalModel;
import repicea.util.ObjectUtility;
import repicea.util.REpiceaLogManager;

public class SimulationStudy {

	static boolean VERBOSE = false;
	
	class WorkerThread extends Thread {
		final SpatialPopulation p;
		final int real;
		final int n;
				
		WorkerThread(int id, SpatialPopulation p, int real, int n) {
			this.setName("Worker " + id);
			this.p = p.clone();
			this.real = real;
			this.n = n;
		}
		
		@Override
		public void run() {
			int moduloValue = VERBOSE ? 1 : 1000;
			int successfulSampling = 0;
			int unsuccessfulSampling = 0;
			while (successfulSampling < real) {
				if (successfulSampling > 0 && successfulSampling%moduloValue == 0) {
					System.out.println(getName() + " is processing realization " + successfulSampling + "(" + unsuccessfulSampling + " failures)");
				}
				Object[] realization = p.samplePopulation(n);
				if (realization != null) {
					synchronized(lock) {
						try {
							writer.addRecord(realization);
							successfulSampling++;
						} catch (IOException e) {
							e.printStackTrace();
							throw new RuntimeException(e.getMessage());
						}
					}
				} else {
					unsuccessfulSampling++;
				}
			}
		}

		
	}


	private Object lock = new Object();
	private final CSVWriter writer;
	
	SimulationStudy(int nbRealizations, int popSize, int sampleSize, int nbThreads) throws IOException, InterruptedException {
		String suffix = "S";
		String path = ObjectUtility.getPackagePath(getClass());
		
		// Set the csv writer
		String filename = path + "simulation_" + suffix + "_N" + popSize * popSize + "_n" + sampleSize + ".csv";
		writer = new CSVWriter(new File(filename), false);
		List<FormatField> fields = new ArrayList<FormatField>();
		fields.add(new CSVField("true_beta0"));
		fields.add(new CSVField("true_beta1"));
		fields.add(new CSVField("beta0"));
		fields.add(new CSVField("var_beta0"));
		fields.add(new CSVField("beta1"));
		fields.add(new CSVField("var_beta1"));
		fields.add(new CSVField("beta0_GLM"));
		fields.add(new CSVField("beta1_GLM"));
		writer.setFields(fields);
		
		SpatialPopulation pop;
		System.out.println("Creating population...");
		// Create a new population
//		Matrix trueBeta = new Matrix(2,1);
//		trueBeta.setValueAt(0, 0, 1.0);
//		trueBeta.setValueAt(1, 0, -0.1);
//		pop = new SpatialPopulation(trueBeta, popSize);
//		String popFilename = path + "population" + suffix + "_" + popSize + ".csv";
//		pop.save(popFilename);
//		XmlSerializer serializer = new XmlSerializer(path + "population" + suffix + "_" + popSize + ".zml");
//		serializer.writeObject(pop);
		
		// Deserialize the population
		XmlDeserializer deserializer = new XmlDeserializer(path + "population" + suffix + "_" + popSize + ".zml");
		pop = (SpatialPopulation) deserializer.readObject();
		System.out.println("Population created.");
		
		
		if (nbThreads > 1) {
			Thread[] threadArray = new Thread[nbThreads];
			int real = (int) ((double) nbRealizations / nbThreads);
			for (int k = 0; k < threadArray.length; k++) {
				Thread t = new WorkerThread(k, pop, real, sampleSize);
				threadArray[k] = t;
				t.start();
			}
			for (Thread t : threadArray) {
				t.join();
			}
		} else {
			new WorkerThread(0, pop, nbRealizations, sampleSize).run();
		}
		writer.close();
	}

	

	public static void main(String[] args) throws Exception {
		AbstractStatisticalModel.LOGGER_NAME = MaximumLikelihoodEstimator.LOGGER_NAME;
		REpiceaLogManager.getLogger(MaximumLikelihoodEstimator.LOGGER_NAME).setLevel(Level.OFF);
		SimulationStudy.VERBOSE = true;
		new SimulationStudy(1000, 200, 2000, 2);
	}
	
}
