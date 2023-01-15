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

import repicea.io.FormatField;
import repicea.io.Saveable;
import repicea.io.javacsv.CSVField;
import repicea.io.javacsv.CSVWriter;
import repicea.stats.QuantileUtility;
import repicea.stats.StatisticalUtility;
import repicea.stats.estimates.MonteCarloEstimate;
import repicea.stats.model.glm.measerr.simulation.DistanceCalculator.SpatialPopulationUnitDistance;


@SuppressWarnings("serial")
class DistanceCalculator extends ArrayList<SpatialPopulationUnitDistance> implements Saveable {
	
	static final class SpatialPopulationUnitDistance implements Comparable<SpatialPopulationUnitDistance> {
		final double distance;
		final SpatialPopulationUnit pu;
		
		SpatialPopulationUnitDistance(double distance, SpatialPopulationUnit pu) {
			this.distance = distance;
			this.pu = pu;
		}

		@Override
		public int compareTo(SpatialPopulationUnitDistance o) {
			if (distance < o.distance) {
				return -1;
			} else if (distance > o.distance) {
				return 1;
			}
			return 0;
		}
		
	}

	final SpatialPopulationUnit thisUnit;
	
	DistanceCalculator(SpatialPopulationUnit pu) {
		thisUnit = pu;
	}
	
	void setDistanceToConspecific(List<SpatialPopulationUnit> pus, boolean isPopulation) {
		clear();
		double p = 0.01;
		List<Double> distances = new ArrayList<Double>();
		for (SpatialPopulationUnit pu : pus) {
			if (pu.isConspecificIn) {
				distances.add(thisUnit.getDistanceFrom(pu));
			}
		}
		if (isPopulation) {
			thisUnit.trueDistanceToConspecific = QuantileUtility.getQuantileFromPopulation(distances, p);
		} else {
			thisUnit.measuredDistanceToConspecific = QuantileUtility.getQuantileEstimateFromSample(distances, p, null);
			MonteCarloEstimate estimate = QuantileUtility.getQuantileEstimateFromSample(distances, p, null, 100);
			thisUnit.variance = estimate.getVariance().getValueAt(0, 0);
		}
	}
	
	List<SpatialPopulationUnit> retrievePopulationUnits() {
		List<SpatialPopulationUnit> newList = new ArrayList<SpatialPopulationUnit>();
		for (SpatialPopulationUnitDistance spud : this) {
			newList.add(spud.pu);
		}
		return newList;
	}

	@Override
	public void save(String filename) throws IOException {
		CSVWriter writer = new CSVWriter(new File(filename), false);
		List<FormatField> fields = new ArrayList<FormatField>();
		fields.add(new CSVField("id"));
		fields.add(new CSVField("distance"));
		fields.add(new CSVField("isConspecificIn"));
		writer.setFields(fields);
		
		Object[] record;
//		int nbOccurrences = 0;
		for (SpatialPopulationUnitDistance spdu : this) {
			record = new Object[3];
			record[0] = spdu.pu.id;
			record[1] = spdu.distance;
			record[2] = spdu.pu.isConspecificIn ? 1 : 0;
			writer.addRecord(record);
//			if (spdu.pu.isConspecificIn) {
//				nbOccurrences++;
//			}
//			if (nbOccurrences >= 30) {
//				break;
//			}
		}
		writer.close();
	}
}
