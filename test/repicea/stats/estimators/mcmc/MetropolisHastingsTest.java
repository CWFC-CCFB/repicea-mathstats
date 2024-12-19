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

import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.FixMethodOrder;
import org.junit.Test;
import org.junit.runners.MethodSorters;

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

@FixMethodOrder(MethodSorters.NAME_ASCENDING)
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
		public int getNbSubjects() {
			return getNumberOfObservations();
		}

		@Override
		public double getLikelihoodOfThisSubject(Matrix parms, int subjectId) {
			return GaussianUtility.getProbabilityDensity(observations.get(subjectId), parms.getValueAt(0, 0), parms.getValueAt(1, 0));
		}

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
		Map<Integer, List<Integer>> hierarchicalMap;
		
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

		private Map<Integer, List<Integer>> getHierarchicalMap() {
			if (hierarchicalMap == null) {
				hierarchicalMap = new HashMap<Integer, List<Integer>>();
				for (int i = 0; i < getNumberOfObservations(); i++) {
					MixedModelObservation o = observations.get(i);
					if (!hierarchicalMap.containsKey(o.subjectID)) {
						hierarchicalMap.put(o.subjectID, new ArrayList<Integer>());
					}
					hierarchicalMap.get(o.subjectID).add(i);
				}
			}
			return hierarchicalMap;
		}
		
		@Override
		public int getNbSubjects() {
			return getHierarchicalMap().size();
		}

		@Override
		public double getLikelihoodOfThisSubject(Matrix parms, int subjectId) {
			List<Integer> observationIndices = getHierarchicalMap().get(subjectId);
			double randomEffect = parms.getValueAt(3 + subjectId, 0);	
			double jointProb = 1d;
			for (int index : observationIndices) {
				MixedModelObservation o = observations.get(index);
				jointProb *= GaussianUtility.getProbabilityDensity(o.observation, parms.getValueAt(0, 0) + randomEffect, parms.getValueAt(1, 0));
			}
			return jointProb;
		}

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
		
		
		MetropolisHastingsCompatibleModelImpl convertToSimpleImplementation() {
			MetropolisHastingsCompatibleModelImpl impl = new MetropolisHastingsCompatibleModelImpl(0);
			for (MixedModelObservation o : observations) {
				impl.observations.add(o.observation);
			}
			return impl;
		}
		
	}

	private static MetropolisHastingsAlgorithm hierarchicalWithoutRandomEffects;
	private static MetropolisHastingsAlgorithm hierarchicalWithRandomEffects;
	
	@Test
	public void test01SimpleMCMC() throws MarshallingException, UnmarshallingException {
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
	public void test02HierarchicalDataWithoutMixedEffectsMCMC() throws MarshallingException, UnmarshallingException {
		MetropolisHastingsCompatibleModel mhcm;
		
		String filename = ObjectUtility.getPackagePath(getClass()) + "refObservationsMixedEffects.zml";
		
//		UNCOMMENT THE NEXT LINES TO UPDATE THE OBSERVATIONS		
//		mhcm = new MetropolisHastingsCompatibleMixedModelImpl(10, 100);
//		XmlSerializer serializer = new XmlSerializer(filename);
//		serializer.writeObject(mhcm);
		
		XmlDeserializer deserializer = new XmlDeserializer(filename);
		MetropolisHastingsCompatibleMixedModelImpl mixedImpl = (MetropolisHastingsCompatibleMixedModelImpl) deserializer.readObject();
		mhcm = mixedImpl.convertToSimpleImplementation();
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
		Assert.assertEquals("Testing variance", 19.5, parameterEstimates.getValueAt(1, 0), 0.1);
		
		DataSet ds = mha.convertMetropolisHastingsSampleToDataSet();
		Assert.assertEquals("Testing field name 0", "LLK", ds.getFieldNames().get(0));
		Assert.assertEquals("Testing field name 1", "Mean", ds.getFieldNames().get(1));
		Assert.assertEquals("Testing field name 2", "Variance", ds.getFieldNames().get(2));
		Assert.assertEquals("Testing posterior parameter sample size", 10000, ds.getNumberOfObservations());
		
		hierarchicalWithoutRandomEffects = mha;
	}

	@Test
	public void test03MixedEffectsMCMC() throws MarshallingException, UnmarshallingException {
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
		Assert.assertEquals("Testing mean", 3.8, parameterEstimates.getValueAt(0, 0), 0.2);
		Assert.assertEquals("Testing variance", 15.9, parameterEstimates.getValueAt(1, 0), 0.2);
		Assert.assertEquals("Testing random effect standard deviation", 2.3, parameterEstimates.getValueAt(2, 0), 0.1);
		
		DataSet ds = mha.convertMetropolisHastingsSampleToDataSet();
		Assert.assertEquals("Testing field name 0", "LLK", ds.getFieldNames().get(0));
		Assert.assertEquals("Testing field name 1", "Mean", ds.getFieldNames().get(1));
		Assert.assertEquals("Testing field name 2", "Variance", ds.getFieldNames().get(2));
		Assert.assertEquals("Testing posterior parameter sample size", 10000, ds.getNumberOfObservations());

		hierarchicalWithRandomEffects = mha;
	}

	@Test
	public void test04ModelComparison() {
		double lpmlWithout = hierarchicalWithoutRandomEffects.getLogPseudomarginalLikelihood();
		double lpmlWith = hierarchicalWithRandomEffects.getLogPseudomarginalLikelihood();
		double pseudoBayesFactor = Math.exp(lpmlWith - lpmlWithout);
		Assert.assertTrue("Testing pseudo Bayes factor greater than 10", pseudoBayesFactor > 10);
	}
	
	@Test
	public void test04NoGridTest() {
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
	public void test05WithGridTest() {
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

	@Test
	public void test06ReleasingFinalSample() {
		DataSet ds1 = hierarchicalWithoutRandomEffects.getParameterEstimatesReport();
		System.out.println(ds1.toString());
		hierarchicalWithoutRandomEffects.releaseFinalSampleSelection();
		DataSet ds2 = hierarchicalWithoutRandomEffects.getParameterEstimatesReport();
		System.out.println(ds2.toString());
		Assert.assertEquals("Comparing the two outputs", ds1.toString(), ds2.toString());
	}

	@AfterClass
    public static void afterClass() {
		hierarchicalWithoutRandomEffects = null;
		hierarchicalWithRandomEffects = null;
    }
}
