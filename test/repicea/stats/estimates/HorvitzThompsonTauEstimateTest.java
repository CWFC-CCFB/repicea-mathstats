/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2025 His Majesty the King in right of Canada
 * Author: Mathieu Fortin, Canadian Forest Service
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
package repicea.stats.estimates;

import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.stats.sampling.PopulationUnit;

public class HorvitzThompsonTauEstimateTest {

	@Test
	public void simpleTotalAndVarianceEstimateTest() {
		List<Double> sample = new ArrayList<Double>();
		sample.add(2d);
		sample.add(4d);
		sample.add(2d);
		sample.add(5d);
		sample.add(7d);
		sample.add(1d);
		sample.add(5d);
		sample.add(4d);
		sample.add(7d);
		PopulationTotalEstimate estimate = new PopulationTotalEstimate(1000);
		Matrix obs;
		for (int i = 0; i < sample.size(); i++) {
			double value = sample.get(i);
			obs = new Matrix(1,1);
			obs.setValueAt(0, 0, value);
			estimate.addObservation(new PopulationUnit(i + "", obs));
		}
		Matrix total = estimate.getMean();
		Assert.assertEquals("Testing the estimate of the total", 4111.11111111111, total.getValueAt(0, 0), 1E-8);
		Matrix totalVariance = estimate.getVariance();
		Assert.assertEquals("Testing the variance of the total", 507734.5679, totalVariance.getValueAt(0, 0), 1E-4);
	}
	
	
	
	
}
