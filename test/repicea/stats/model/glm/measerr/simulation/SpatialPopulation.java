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

import java.util.List;

import repicea.math.Matrix;
import repicea.stats.StatisticalUtility;
import repicea.stats.data.DataSet;
import repicea.stats.estimates.Estimate;
import repicea.stats.model.glm.GeneralizedLinearModel;
import repicea.stats.model.glm.Family.GLMDistribution;
import repicea.stats.model.glm.LinkFunction.Type;
import repicea.stats.model.glm.measerr.SIMEXModel;
import repicea.stats.sampling.SamplingUtility;
import repicea.util.ObjectUtility;

class SpatialPopulation extends AbstractPopulation<SpatialPopulationUnit> {

	static String PATH = ObjectUtility.getPackagePath(SpatialPopulation.class);
	static int realID = 0;
	
	SpatialPopulation(Matrix trueBeta, int size) {
		this(trueBeta, size, false);
	}
	
	private SpatialPopulation(Matrix trueBeta, int size, boolean isClone) {
		super(trueBeta, size);
		if (!isClone) {
			int id = 0;
			for (int x = 1; x <= size; x++) {
				for (int y = 1; y <= size; y++) {
					if (y <= size * .5) {
						populationUnits.add(new SpatialPopulationUnit(id++, x, y, StatisticalUtility.getRandom().nextDouble() < 0.6));
					} else {
						populationUnits.add(new SpatialPopulationUnit(id++, x, y, false));
					}
				}
			}
			determineDistanceToConspecific(populationUnits, true);	// determine the nearest conspecific throughout the population units
			generateRealizations();
		}
	}

	@Override
	public SpatialPopulation clone() {
		SpatialPopulation pop = new SpatialPopulation(trueBeta, size, true);
		for (SpatialPopulationUnit pu : populationUnits) {
			pop.populationUnits.add(pu.clone());
		}
		return pop;
	}
	
	static void determineDistanceToConspecific(List<SpatialPopulationUnit> popUnits, boolean isPopulation) {
		for (int i = 0; i < popUnits.size(); i++) {
			SpatialPopulationUnit pu_i = popUnits.get(i);
			pu_i.dc.setDistanceToConspecific(popUnits, isPopulation);
		}
	}

	@Override
	public Object[] samplePopulation(int n) {
		Object[] record = new Object[6];
		int nbPositive = generateRealizations();
//		System.out.println("Nb positive cases = " + nbPositive);
		List<SpatialPopulationUnit> sample = SamplingUtility.getSample(populationUnits, n);
		int nbPositiveInSample = 0;
		for (SpatialPopulationUnit spu : sample) {
			nbPositiveInSample += spu.getY();
		}
		determineDistanceToConspecific(sample, false); // false is the sample
		
		DataSet ds = new DataSet(sample.get(0).getFieldname(false));
		for (PopulationUnit pu : sample) {
			ds.addObservation(pu.asObservation(false));
//			if (!((SpatialPopulationUnit) pu).isConspecificIn) {
//				try {
//					((SpatialPopulationUnit) pu).dc.save(PATH + "exampleDistanceSample.csv");
//				} catch (Exception e) {}
//			}
		}
		ds.indexFieldType();
		try {
			if (realID < 10)
				ds.save(PATH + "sample" + realID++ + ".csv");
		} catch(Exception e) {}
		
		GeneralizedLinearModel pre_glm = new GeneralizedLinearModel(ds,  GLMDistribution.Bernoulli, Type.CLogLog, "y ~ distanceToConspecific");
		pre_glm.doEstimation();
//		pre_glm.getSummary();
//		GLMWithMeasurementError glm = new GLMWithMeasurementError(ds, 
//				"y ~ distanceToConspecific", 
//				pre_glm.getParameters(),
//				new GLMNormalClassicalMeasErrorDefinition("distanceToConspecific", "variance", .1));
//		((MaximumLikelihoodEstimator) glm.getEstimator()).setLineSearchMethod(LineSearchMethod.HALF_STEP);
		SIMEXModel glm = new SIMEXModel(pre_glm, "distanceToConspecific", "variance");
		glm.doEstimation();
//		glm.getSummary();
		
		if (glm.getEstimator().isConvergenceAchieved()) {
			record = new Object[8];
			Estimate<?> est = glm.getEstimator().getParameterEstimates();
			record[0] = trueBeta.getValueAt(0, 0);
			record[1] = trueBeta.getValueAt(1, 0);
			record[2] = est.getMean().getValueAt(0, 0);
			record[3] = est.getVariance().getValueAt(0, 0);
			record[4] = est.getMean().getValueAt(1, 0);
			record[5] = est.getVariance().getValueAt(1, 1);
			record[6] = pre_glm.getParameters().getValueAt(0, 0);
			record[7] = pre_glm.getParameters().getValueAt(1, 0);
			return record;
		} else {
			return null;
		}

	}

}
