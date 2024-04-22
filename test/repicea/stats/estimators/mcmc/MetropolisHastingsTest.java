/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2024 His Majesty the King in Right of Canada
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
package repicea.stats.estimators.mcmc;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.math.utility.GaussianUtility;
import repicea.serial.MarshallingException;
import repicea.serial.UnmarshallingException;
import repicea.serial.xml.XmlDeserializer;
import repicea.stats.StatisticalUtility;
import repicea.stats.data.DataSet;
import repicea.stats.distributions.ContinuousDistribution;
import repicea.stats.distributions.GaussianDistribution;
import repicea.stats.distributions.UniformDistribution;
import repicea.util.ObjectUtility;

public class MetropolisHastingsTest {

	static class MetropolisHastingsCompatibleModelImpl implements MetropolisHastingsCompatibleModel {

		static final double TrueMean = 3;
		static final double TrueSTD = 4;
		
		final List<Double> observations;
		
		MetropolisHastingsCompatibleModelImpl(int size) {
			observations = new ArrayList<Double>();
			for (int i = 0; i < size; i++) {
				observations.add(StatisticalUtility.getRandom().nextGaussian() * TrueSTD + TrueMean);
			}
		}
		
		
		@Override
		public boolean isInterceptModel() {return false;}

		@Override
		public List<String> getEffectList() {
			List<String> effects = new ArrayList<String>();
			effects.add("Mean");
			effects.add("Variance");
			return effects;
		}

		@Override
		public int getNumberOfObservations() {
			return observations.size();
		}

		@Override
		public List<String> getOtherParameterNames() {
			return new ArrayList<String>();
		}

		@Override
		public double getLogLikelihood(Matrix parms) {
			double llk = 0;
			for (int i = 0; i < getNbSubjects(); i++) {
				llk += Math.log(GaussianUtility.getProbabilityDensity(observations.get(i), parms.getValueAt(0, 0), parms.getValueAt(1, 0)));
			}
			return llk;
		}

		@Override
		public int getNbSubjects() {
			return getNumberOfObservations();
		}

		// useless
		@Override
		public double getLikelihoodOfThisSubject(Matrix parms, int subjectId) {return 0d;}

		@Override
		public GaussianDistribution getStartingParmEst(double coefVar) {
			Matrix parmEst = new Matrix(2,1);
			parmEst.setValueAt(0, 0, 1d);
			parmEst.setValueAt(1, 0, 10d);

			Matrix varianceDiag = new Matrix(parmEst.m_iRows,1);
			for (int i = 0; i < varianceDiag.m_iRows; i++) {
				varianceDiag.setValueAt(i, 0, Math.pow(parmEst.getValueAt(i, 0) * coefVar, 2d));
			}
			
			GaussianDistribution gd = new GaussianDistribution(parmEst, varianceDiag.matrixDiagonal());

			return gd;
		}

		@Override
		public void setPriorDistributions(MetropolisHastingsPriorHandler handler) {
			handler.addFixedEffectDistribution(new UniformDistribution(-10, 10), 0);
			handler.addFixedEffectDistribution(new UniformDistribution(0, 25), 1);
		}	
		
	}

	static class MixedModelObservation {
		final int subjectID;
		final double observation;
		MixedModelObservation(int subjectID, double observation) {
			this.subjectID = subjectID;
			this.observation = observation;
		}
	}
	
	static class MetropolisHastingsCompatibleMixedModelImpl implements MetropolisHastingsCompatibleModel {

		static final double TrueMean = 3;
		static final double TrueVariance = 16;
		static final double TrueSTD_u = 2;
		final List<MixedModelObservation> observations;
		final Map<Integer, Double> randomEffects;
		
		MetropolisHastingsCompatibleMixedModelImpl(int nbSubjects, int nbObsPerSubject) {
			double trueSTD = Math.sqrt(TrueVariance);
			randomEffects = new HashMap<Integer, Double>();
			observations = new ArrayList<MixedModelObservation>();
			for (int i = 0; i < nbSubjects; i++) {
				double randomEffect = StatisticalUtility.getRandom().nextGaussian() * TrueSTD_u;
				randomEffects.put(i, randomEffect);
				for (int j = 0; j < nbObsPerSubject; j++) {
					double obs = StatisticalUtility.getRandom().nextGaussian() * trueSTD + TrueMean + randomEffect;
					observations.add(new MixedModelObservation(i, obs));
				}
			}
		}
		
		
		@Override
		public boolean isInterceptModel() {return false;}

		@Override
		public List<String> getEffectList() {
			List<String> effects = new ArrayList<String>();
			effects.add("Mean");
			effects.add("Variance");
			return effects;
		}

		@Override
		public int getNumberOfObservations() {
			return observations.size();
		}

		@Override
		public List<String> getOtherParameterNames() {
			List<String> otherEffects = new ArrayList<String>();
			otherEffects.add("sigma_u");
			for (int i = 0; i < randomEffects.size(); i++) {
				otherEffects.add("u_" + i);
			}
			return otherEffects;
		}

		@Override
		public double getLogLikelihood(Matrix parms) {
			double llk = 0;
			for (MixedModelObservation o : observations) {
				double randomEffect = parms.getValueAt(3 + o.subjectID, 0);	
				llk += Math.log(GaussianUtility.getProbabilityDensity(o.observation, parms.getValueAt(0, 0) + randomEffect, parms.getValueAt(1, 0)));
			}
			return llk;
		}

		@Override
		public int getNbSubjects() {
			return getNumberOfObservations();
		}

		// useless
		@Override
		public double getLikelihoodOfThisSubject(Matrix parms, int subjectId) {return 0d;}

		@Override
		public GaussianDistribution getStartingParmEst(double coefVar) {
			Matrix parmEst = new Matrix(3 + randomEffects.size(),1);
			parmEst.setValueAt(0, 0, 1d);
			parmEst.setValueAt(1, 0, 10d);
			parmEst.setValueAt(2, 0, 5d);

			for (int i = 0; i < randomEffects.size(); i++) {
				parmEst.setValueAt(3 + i, 0, 0);
			}

			Matrix varianceDiag = new Matrix(parmEst.m_iRows,1);
			for (int i = 0; i < varianceDiag.m_iRows; i++) {
				double varianceSampler = i < 3 ?
						Math.pow(parmEst.getValueAt(i, 0) * coefVar, 2d) :
							Math.pow(parmEst.getValueAt(2, 0) * coefVar, 2d);

				varianceDiag.setValueAt(i, 0, varianceSampler);
			}
			
			GaussianDistribution gd = new GaussianDistribution(parmEst, varianceDiag.matrixDiagonal());

			return gd;
		}

		@Override
		public void setPriorDistributions(MetropolisHastingsPriorHandler handler) {
			handler.addFixedEffectDistribution(new UniformDistribution(-10, 10), 0);
			handler.addFixedEffectDistribution(new UniformDistribution(0, 25), 1);
			ContinuousDistribution stdPrior = new UniformDistribution(0, 10);
			handler.addFixedEffectDistribution(stdPrior, 2);
			for (int i = 0; i < randomEffects.size(); i++) {
				handler.addRandomEffectStandardDeviation(new GaussianDistribution(0, 1), stdPrior, 3 + i);
			}
		}	
		
	}

	@Test
	public void simpleMCMC() throws MarshallingException, UnmarshallingException {
		MetropolisHastingsCompatibleModel mhcm;
		
		String filename = ObjectUtility.getPackagePath(getClass()) + "refObservations.zml";
		
//		UNCOMMENT THE NEXT LINES TO UPDATE THE OBSERVATIONS		
//		mhcm = new MetropolisHastingsCompatibleModelImpl(100);
//		XmlSerializer serializer = new XmlSerializer(filename);
//		serializer.writeObject(mhcm);
		
		XmlDeserializer deserializer = new XmlDeserializer(filename);
		mhcm = (MetropolisHastingsCompatibleModel) deserializer.readObject();
		
		MetropolisHastingsAlgorithm mha = new MetropolisHastingsAlgorithm(mhcm, "mh", "MH");
		System.out.println(mha.getSimulationParameters().toString());
//		Level l = Level.FINE;
//		ConsoleHandler ch = new ConsoleHandler();
//		ch.setLevel(l);
//		REpiceaLogManager.getLogger("mh").setLevel(Level.FINE);
//		REpiceaLogManager.getLogger("mh").addHandler(ch);

		mha.doEstimation();
		System.out.println(mha.getReport());
		Matrix parameterEstimates = mha.getFinalParameterEstimates();
		Assert.assertEquals("Testing mean", 2.49, parameterEstimates.getValueAt(0, 0), 0.1);
		Assert.assertEquals("Testing variance", 15.5, parameterEstimates.getValueAt(1, 0), 0.1);
		
		DataSet ds = mha.convertMetropolisHastingsSampleToDataSet();
		Assert.assertEquals("Testing field name 0", "LLK", ds.getFieldNames().get(0));
		Assert.assertEquals("Testing field name 1", "Mean", ds.getFieldNames().get(1));
		Assert.assertEquals("Testing field name 2", "Variance", ds.getFieldNames().get(2));
		Assert.assertEquals("Testing posterior parameter sample size", 10000, ds.getNumberOfObservations());
	}

	
	@Test
	public void mixedEffectsMCMC() throws MarshallingException, UnmarshallingException {
		MetropolisHastingsCompatibleModel mhcm;
		
		String filename = ObjectUtility.getPackagePath(getClass()) + "refObservationsMixedEffects.zml";
		
//		UNCOMMENT THE NEXT LINES TO UPDATE THE OBSERVATIONS		
//		mhcm = new MetropolisHastingsCompatibleMixedModelImpl(10, 100);
//		XmlSerializer serializer = new XmlSerializer(filename);
//		serializer.writeObject(mhcm);
		
		XmlDeserializer deserializer = new XmlDeserializer(filename);
		mhcm = (MetropolisHastingsCompatibleMixedModelImpl) deserializer.readObject();
		
		MetropolisHastingsAlgorithm mha = new MetropolisHastingsAlgorithm(mhcm, "mh", "MH");
		System.out.println(mha.getSimulationParameters().toString());
//		Level l = Level.FINE;
//		ConsoleHandler ch = new ConsoleHandler();
//		ch.setLevel(l);
//		REpiceaLogManager.getLogger("mh").setLevel(Level.FINE);
//		REpiceaLogManager.getLogger("mh").addHandler(ch);

		mha.doEstimation();
		System.out.println(mha.getReport());
		Matrix parameterEstimates = mha.getFinalParameterEstimates();
		Assert.assertEquals("Testing mean", 3.8, parameterEstimates.getValueAt(0, 0), 0.1);
		Assert.assertEquals("Testing variance", 15.9, parameterEstimates.getValueAt(1, 0), 0.1);
		Assert.assertEquals("Testing random effect standard deviation", 2.3, parameterEstimates.getValueAt(2, 0), 0.1);
		
		DataSet ds = mha.convertMetropolisHastingsSampleToDataSet();
		Assert.assertEquals("Testing field name 0", "LLK", ds.getFieldNames().get(0));
		Assert.assertEquals("Testing field name 1", "Mean", ds.getFieldNames().get(1));
		Assert.assertEquals("Testing field name 2", "Variance", ds.getFieldNames().get(2));
		Assert.assertEquals("Testing posterior parameter sample size", 10000, ds.getNumberOfObservations());
	}

	
	@Test
	public void noGridTest() {
		MetropolisHastingsCompatibleModel model = new MetropolisHastingsCompatibleModelImpl(100);
		MetropolisHastingsAlgorithm mha = new MetropolisHastingsAlgorithm(model, "mh", "MH");
		model.setPriorDistributions(mha.getPriorHandler());
		mha.simParms.nbInitialGrid = 0;
		double coefVar = 0.01;
		GaussianDistribution samplingDist = model.getStartingParmEst(coefVar);
		MetropolisHastingsSample firstSet = mha.getFirstSetOfParameters(samplingDist);
		GaussianDistribution refSamplingDist = model.getStartingParmEst(coefVar);
		Assert.assertEquals("Testing first parameter", 
				refSamplingDist.getMean().getValueAt(0,0),
				firstSet.parms.getValueAt(0, 0), 1E-8);
		Assert.assertEquals("Testing second parameter", 
				refSamplingDist.getMean().getValueAt(1,0),
				firstSet.parms.getValueAt(1, 0), 1E-8);
	}

	@Test
	public void withGridTest() {
		MetropolisHastingsCompatibleModel model = new MetropolisHastingsCompatibleModelImpl(100);
		MetropolisHastingsAlgorithm mha = new MetropolisHastingsAlgorithm(model, "mh", "MH");
		model.setPriorDistributions(mha.getPriorHandler());
		mha.simParms.nbInitialGrid = 1000;
		double coefVar = 0.01;
		GaussianDistribution samplingDist = model.getStartingParmEst(coefVar);
		MetropolisHastingsSample firstSet = mha.getFirstSetOfParameters(samplingDist);
		GaussianDistribution refSamplingDist = model.getStartingParmEst(coefVar);
		Assert.assertTrue("Testing first parameter", Math.abs(
				refSamplingDist.getMean().getValueAt(0,0) - firstSet.parms.getValueAt(0, 0)) > 1E-8);
		Assert.assertTrue("Testing secondirst parameter", Math.abs(
				refSamplingDist.getMean().getValueAt(0,0) - firstSet.parms.getValueAt(0, 0)) > 1E-8);
	}

}
