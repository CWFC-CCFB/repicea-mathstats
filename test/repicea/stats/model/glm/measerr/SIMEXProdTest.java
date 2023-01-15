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
package repicea.stats.model.glm.measerr;

import java.util.logging.ConsoleHandler;
import java.util.logging.Level;

import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import repicea.math.optimizer.NewtonRaphsonOptimizer;
import repicea.stats.data.DataSet;
import repicea.stats.estimators.MaximumLikelihoodEstimator;
import repicea.stats.model.glm.GeneralizedLinearModel;
import repicea.stats.model.glm.Family.GLMDistribution;
import repicea.stats.model.glm.LinkFunction.Type;
import repicea.util.ObjectUtility;
import repicea.util.REpiceaLogManager;

public class SIMEXProdTest {

	public static void main(String[] args) throws Exception {
		Level l = Level.OFF;
		NewtonRaphsonOptimizer.LOGGER_NAME = MaximumLikelihoodEstimator.LOGGER_NAME;
		ConsoleHandler ch = new ConsoleHandler();
		ch.setLevel(l);
		REpiceaLogManager.getLogger(MaximumLikelihoodEstimator.LOGGER_NAME).setLevel(l);
		REpiceaLogManager.getLogger(MaximumLikelihoodEstimator.LOGGER_NAME).addHandler(ch);		
 		String filename = ObjectUtility.getPackagePath(GLModelWithMeasErrorTest.class).concat("OccurrencePartDataSet_ERS.csv");
		DataSet dataSet = new DataSet(filename, true);
		String formula = "occurred ~ lnDt + DD + logDD + TotalPrcp + logPrcp + LowestTmin + dummyDrainage4hydrique + G_F + lnG_F + lnG_R + occIndex10km + sqr(occIndex10km) + speciesThere + lnG_SpGr + timeSince1970";
//		String formula = "occurred ~ lnDt + DD + logDD + TotalPrcp + logPrcp + LowestTmin + dummyDrainage4hydrique + G_F + lnG_F + lnG_R + occIndex10km + speciesThere + lnG_SpGr + timeSince1970";
		GeneralizedLinearModel glm = new GeneralizedLinearModel(dataSet, GLMDistribution.Bernoulli, Type.CLogLog, formula);
		glm.doEstimation();
		System.out.print(glm.getSummary());
		SIMEXModel s = new SIMEXModel(glm, "occIndex10km", "occIndex10kmVar");
		s.setNumberOfBootstrapRealizations(100);
		s.setNbThreads(4);
		s.doEstimation();
		System.out.print(s.getSummary());
//		Assert.assertTrue("Checking if successfully extrapolated", s.getEstimator().isConvergenceAchieved());
//		Estimate<?> estimate = s.getEstimator().getParameterEstimates();
//		Assert.assertEquals("Checking parm estimate 1", -0.1687, estimate.getMean().getValueAt(0, 0), 4E-3);
//		Assert.assertEquals("Checking standard error parm estimate 1", 0.067178, Math.sqrt(estimate.getVariance().getValueAt(0, 0)), 4E-4);
//		Assert.assertEquals("Checking parm 2", -0.05505, estimate.getMean().getValueAt(1, 0), 4E-4);
//		Assert.assertEquals("Checking standard error parm estimate 2", 0.004364, Math.sqrt(estimate.getVariance().getValueAt(1, 1)), 4E-5);
	}

}
