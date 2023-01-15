/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2016 Mathieu Fortin for Rouge-Epicea
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
package repicea.simulation;

import java.io.Serializable;
import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CopyOnWriteArrayList;

import repicea.math.Matrix;
import repicea.simulation.REpiceaPredictorEvent.ModelBasedSimulatorEventProperty;
import repicea.stats.distributions.GaussianErrorTerm;
import repicea.stats.distributions.GaussianErrorTermList;
import repicea.stats.distributions.GaussianErrorTermList.IndexableErrorTerm;
import repicea.stats.distributions.StandardGaussianDistribution;
import repicea.stats.estimates.Estimate;
import repicea.stats.estimates.GaussianErrorTermEstimate;

/**
 * This class is the basic class for all models that are designed for predictions.
 * It implements the methods for stochastic simulations.
 * @author Mathieu Fortin - October 2011
 */
@SuppressWarnings("serial")
public abstract class REpiceaPredictor extends SensitivityAnalysisParameter<ModelParameterEstimates> {

	protected static final List<Integer> DefaultZeroIndex = new ArrayList<Integer>();
	static {
		DefaultZeroIndex.add(0);
	}
	
	/**
	 * This class creates a fake subject for interval random effects nested in the plots. 
	 * @author Mathieu Fortin - November 2016
	 */
	protected static class IntervalNestedInPlotDefinition implements MonteCarloSimulationCompliantObject, Serializable {

		private final int monteCarloRealizationID;
		private final String subjectID;
		
		protected IntervalNestedInPlotDefinition(MonteCarloSimulationCompliantObject stand, int date) {
			monteCarloRealizationID = stand.getMonteCarloRealizationId();
			subjectID = getSubjectID(stand, date);
		}
		
		
		@Override
		public String getSubjectId() {
			return subjectID;
		}

		@Override
		public HierarchicalLevel getHierarchicalLevel() {
			return HierarchicalLevel.INTERVAL_NESTED_IN_PLOT;
		}


		@Override
		public int getMonteCarloRealizationId() {
			return monteCarloRealizationID;
		}
		
		private static String getSubjectID(MonteCarloSimulationCompliantObject stand, int date) {
			return stand.getSubjectId() + "_" + date;
		}
	}

	/**
	 * This class creates a fake subject for cruise line random effects.
	 * @author Mathieu Fortin - April 2017
	 */
	protected static class CruiseLine implements MonteCarloSimulationCompliantObject, Serializable {
		private final String subjectID;
		private final int monteCarloRealizationID;
		
		protected CruiseLine(String subjectID, MonteCarloSimulationCompliantObject subject) {
			this.subjectID = subjectID;
			monteCarloRealizationID = subject.getMonteCarloRealizationId();
		}
		
		@Override
		public String getSubjectId() {return subjectID;}

		@Override
		public HierarchicalLevel getHierarchicalLevel() {return HierarchicalLevel.CRUISE_LINE;}

		@Override
		public int getMonteCarloRealizationId() {return monteCarloRealizationID;}
	}

	
	public static enum ErrorTermGroup {Default}
	
	
	protected final CopyOnWriteArrayList<REpiceaPredictorListener> listeners;
	
//	private boolean areBlupsEstimated;
	
	private final Map<String, CruiseLine> cruiseLineMap;
	private final Map<String, IntervalNestedInPlotDefinition> intervalLists;


	// set by the constructor
	protected final boolean isRandomEffectsVariabilityEnabled;
	protected final boolean isResidualVariabilityEnabled;
		
	protected Matrix oXVector;

	final Map<String, Estimate<? extends StandardGaussianDistribution>> defaultRandomEffects;
	final Map<String, Map<String, Estimate<? extends StandardGaussianDistribution>>> blupsRandomEffects; // key1: hierarchical level, key2: subject id
	final Map<String, List<String>> subjectTestedForBlups; // key: hierarchical level
	
	private final Map<String, Map<String, Matrix>> simulatedRandomEffects;	// refers to the subject + realization ids

	private final Map<Enum<?>, GaussianErrorTermEstimate> defaultResidualError;
	final Map<String, GaussianErrorTermList> simulatedResidualError;		// refers to the subject + realization ids
	
//	protected REpiceaRandom random = new REpiceaRandom();
	
	
	/**
	 * General constructor for all combinations of uncertainty sources.
	 * @param isParametersVariabilityEnabled a boolean that enables the variability at the parameter level
	 * @param isRandomEffectsVariabilityEnabled a boolean that enables the variability at the random effect level
	 * @param isResidualVariabilityEnabled a boolean that enables the variability at the tree level
	 */
	protected REpiceaPredictor(boolean isParametersVariabilityEnabled,
			boolean isRandomEffectsVariabilityEnabled,
			boolean isResidualVariabilityEnabled) {
		super(isParametersVariabilityEnabled);
		this.isRandomEffectsVariabilityEnabled = isRandomEffectsVariabilityEnabled;
		this.isResidualVariabilityEnabled = isResidualVariabilityEnabled;
		
		defaultRandomEffects = new HashMap<String, Estimate<? extends StandardGaussianDistribution>>();
		blupsRandomEffects = new HashMap<String, Map<String, Estimate<? extends StandardGaussianDistribution>>>();
		subjectTestedForBlups = new HashMap<String, List<String>>();		
		
		simulatedRandomEffects = new HashMap<String, Map<String, Matrix>>();
		simulatedResidualError = new HashMap<String, GaussianErrorTermList>();
		
		intervalLists = new HashMap<String, IntervalNestedInPlotDefinition>();
		cruiseLineMap = new HashMap<String, CruiseLine>();

		defaultResidualError = new HashMap<Enum<?>, GaussianErrorTermEstimate>();
		
		listeners = new CopyOnWriteArrayList<REpiceaPredictorListener>();
	}
	
	/**
	 * This method reads all the parameters in .csv files and stores the estimates into members defaultBeta, defaultResidualError,
	 * and defaultRandomEffects.
	 */
	protected abstract void init();
	
	protected Map<String, Estimate<? extends StandardGaussianDistribution>> getDefaultRandomEffects() {
		return defaultRandomEffects;
	}

	@Override
	protected void setParameterEstimates(ModelParameterEstimates gaussianEstimate) {
//		super.setParameterEstimates(new ModelParameterEstimates(gaussianEstimate, this));
		super.setParameterEstimates(gaussianEstimate);
		fireModelBasedSimulatorEvent(new REpiceaPredictorEvent(ModelBasedSimulatorEventProperty.DEFAULT_BETA_JUST_SET, null, getParameterEstimates(), this));
	}
	
//	@Override
//	protected ModelParameterEstimates getParameterEstimates() {
//		return (ModelParameterEstimates) super.getParameterEstimates();
//	}

	protected void setDefaultRandomEffects(HierarchicalLevel level, Estimate<? extends StandardGaussianDistribution> newEstimate) {
		Estimate<? extends StandardGaussianDistribution> formerEstimate = defaultRandomEffects.get(level.getName());
		defaultRandomEffects.put(level.getName(), newEstimate);
		fireModelBasedSimulatorEvent(new REpiceaPredictorEvent(ModelBasedSimulatorEventProperty.DEFAULT_RANDOM_EFFECT_AT_THIS_LEVEL_JUST_SET, null, new Object[]{level, formerEstimate, newEstimate}, this));
	}
	
	protected Estimate<? extends StandardGaussianDistribution> getDefaultRandomEffects(HierarchicalLevel level) {return defaultRandomEffects.get(level.getName());}
	
	protected void setDefaultResidualError(Enum<?> enumVar, GaussianErrorTermEstimate estimate) {
		defaultResidualError.put(enumVar, estimate);
		fireModelBasedSimulatorEvent(new REpiceaPredictorEvent(ModelBasedSimulatorEventProperty.DEFAULT_RESIDUAL_ERROR_JUST_SET, null, new Object[]{enumVar, estimate}, this));
	}
	
	protected GaussianErrorTermEstimate getDefaultResidualError(Enum<?> enumVar) {
		return defaultResidualError.get(enumVar);
	}
	
	
//	/**
//	 * This method generates a stand-specific vector of model parameters using matrix Omega.
//	 * @param subject a MonteCarloSimulationCompliantObject object
//	 */
//	private void setSpecificParametersDeviateForThisRealization(MonteCarloSimulationCompliantObject subject) {
////		getParameterEstimates().simulateBlups(subject);
//	}

	/**
	 * This method checks if the interval definition is available for the stand at that date. If it is, it returns the
	 * instance. Otherwise, it creates a new interval definition.
	 * @param stand A MonteCarloSimulationCompliantObject that designates the stand
	 * @param date an Integer
	 * @return an IntervalDefinition instance
	 */
	protected synchronized IntervalNestedInPlotDefinition getIntervalNestedInPlotDefinition(MonteCarloSimulationCompliantObject stand, int date) {
		String subjectID = IntervalNestedInPlotDefinition.getSubjectID(stand, date);
		String intervalID = getSubjectPlusMonteCarloSpecificId(subjectID, stand.getMonteCarloRealizationId());
		IntervalNestedInPlotDefinition intDef = intervalLists.get(intervalID);
		if (intDef == null) {
			intDef = new IntervalNestedInPlotDefinition(stand, date);
			intervalLists.put(getSubjectPlusMonteCarloSpecificId(intDef), intDef);
		}
		return intDef;
	}

	/**
	 * This method checks if a cruise line exists for this plot
	 * @param cruiseLineID the id of the cruise line
	 * @param stand a MonteCarloSimulationCompliantObject instance
	 * @return a CruiseLine instance
	 */
	protected synchronized CruiseLine getCruiseLineForThisSubject(String cruiseLineID, MonteCarloSimulationCompliantObject stand) {
		String cruiseLineIDPlusMCRealization = cruiseLineID.concat("_") + stand.getMonteCarloRealizationId();
		if (!cruiseLineMap.containsKey(cruiseLineIDPlusMCRealization)) {
			cruiseLineMap.put(cruiseLineIDPlusMCRealization, new CruiseLine(cruiseLineID, stand));
		}
		return cruiseLineMap.get(cruiseLineIDPlusMCRealization);
	}
	
	/**
	 * This method calls the setSpecificParametersDeviateForThisRealization method if the parameter variability is enabled and returns 
	 * a realization-specific simulated vector of model parameters. Otherwise it returns a default vector (beta). Note that the simulated
	 * parameters are related to the Monte Carlo realization. For instance, all subject in a given Monte Carlo realization will have the
	 * same simulation parameters. 
	 * @param subject a subject that implements the MonteCarloSimulationCompliantObject interface
	 * @return a vector of parameters
	 */
	@Override
	protected final synchronized Matrix getParametersForThisRealization(MonteCarloSimulationCompliantObject subject) {
		if (isParametersVariabilityEnabled) {
			if (!simulatedParameters.containsKey(subject.getMonteCarloRealizationId())) {		// the simulated parameters remain constant within the same Monte Carlo iteration
				simulatedParameters.put(subject.getMonteCarloRealizationId(), getParameterEstimates().getRandomDeviate());
				//			setSpecificParametersDeviateForThisRealization(subject);
			}
			return simulatedParameters.get(subject.getMonteCarloRealizationId());
		} else {
			return getParameterEstimates().getMean();
		}
	}

	
	/**
	 * This method generates a subject-specific random effects vector using matrix G.
	 * @param subject a MonteCarloSimulationCompliantObject instance
	 */
	private void setSpecificRandomEffectsForThisSubject(MonteCarloSimulationCompliantObject subject) {
		HierarchicalLevel subjectLevel = subject.getHierarchicalLevel();
		
		Matrix randomDeviates;
		Estimate<? extends StandardGaussianDistribution> originalRandomEffects;
		if (doBlupsExistForThisSubject(subject)) {
			simulateDeviatesForRandomEffectsOfThisSubject(subject, getBlupsForThisSubject(subject));
		} else {
			randomDeviates = simulateDeviatesForRandomEffectsOfThisSubject(subject, defaultRandomEffects.get(subjectLevel.getName()));
			originalRandomEffects = getDefaultRandomEffects(subjectLevel);
			fireRandomEffectDeviateGeneratedEvent(subject, originalRandomEffects, randomDeviates);
		}
	}
	
	protected void fireRandomEffectDeviateGeneratedEvent(MonteCarloSimulationCompliantObject subject,
			Estimate<? extends StandardGaussianDistribution> originalRandomEffects,
			Matrix randomDeviates) {
		REpiceaPredictorEvent event = new REpiceaPredictorEvent(ModelBasedSimulatorEventProperty.RANDOM_EFFECT_DEVIATE_JUST_GENERATED, 
				null, 
				new Object[]{subject, originalRandomEffects, randomDeviates.getDeepClone()},
				this);
		fireModelBasedSimulatorEvent(event);
	}
	
	
	/**
	 * This method simulates random deviates from an estimate and stores them in the simulatedRandomEffects
	 * member.
	 * @param subject a MonteCarloSimulationCompliantObject instance
	 * @param randomEffectsEstimate the estimate from which the random deviates are generated
	 * @return the random deviates as a Matrix instance (a copy of it)
	 */
	protected Matrix simulateDeviatesForRandomEffectsOfThisSubject(MonteCarloSimulationCompliantObject subject, Estimate<?> randomEffectsEstimate) {
		Matrix randomDeviates = randomEffectsEstimate.getRandomDeviate();
		setDeviatesForRandomEffectsOfThisSubject(subject, randomDeviates);
		return randomDeviates.getDeepClone();
	}

	protected final synchronized void setDeviatesForRandomEffectsOfThisSubject(MonteCarloSimulationCompliantObject subject, Matrix randomDeviates) {
		HierarchicalLevel subjectLevel = subject.getHierarchicalLevel();
		if (!simulatedRandomEffects.containsKey(subjectLevel.getName())) {
			simulatedRandomEffects.put(subjectLevel.getName(), new HashMap<String, Matrix>());
		}
		Map<String, Matrix> randomEffectsMap = simulatedRandomEffects.get(subjectLevel.getName());
		randomEffectsMap.put(getSubjectPlusMonteCarloSpecificId(subject), randomDeviates);
	}
	
	
	protected final String getSubjectPlusMonteCarloSpecificId(MonteCarloSimulationCompliantObject object) {
		return getSubjectPlusMonteCarloSpecificId(object.getSubjectId(), object.getMonteCarloRealizationId());
	}

	protected static String getSubjectPlusMonteCarloSpecificId(String subjectID, int monteCarloRealizationID) {
		return subjectID + "_" + monteCarloRealizationID;
	}
	
	/**
	 * This method calls the setSpecificPlotRandomEffectsForThisStand method if the random effects variability is enabled and returns 
	 * a stand-specific simulated vector of random effects. Otherwise it returns a default vector (all elements set to 0).
	 * @param subject a MonteCarloSimulationCompliantObject object
	 * @return a Matrix object
	 */
	protected final synchronized Matrix getRandomEffectsForThisSubject(MonteCarloSimulationCompliantObject subject) {
		HierarchicalLevel subjectLevel = subject.getHierarchicalLevel();
		if (isRandomEffectsVariabilityEnabled) {
			if (!doRandomDeviatesExistForThisSubject(subject)) {
				setSpecificRandomEffectsForThisSubject(subject);
			}
			return simulatedRandomEffects.get(subjectLevel.getName()).get(getSubjectPlusMonteCarloSpecificId(subject));
		} else {
			Estimate<? extends StandardGaussianDistribution> blups = getBlupsForThisSubject(subject);
			if (blups != null) {
				return blups.getMean();
			} else {
				return defaultRandomEffects.get(subjectLevel.getName()).getMean();
			}
		} 
	}
	
	protected final boolean doRandomDeviatesExistForThisSubject(MonteCarloSimulationCompliantObject subject) {
		HierarchicalLevel subjectLevel = subject.getHierarchicalLevel();
		return simulatedRandomEffects.get(subjectLevel.getName()) != null && simulatedRandomEffects.get(subjectLevel.getName()).containsKey(getSubjectPlusMonteCarloSpecificId(subject)); 
	}
	

//	/**
//	 * This method clears all simulated deviates for the parameter estimates, the random effects and the residual errors. 
//	 * IMPORTANT: it does not reset the blups of the random effects. 
//	 */
// 	FIXED - THIS METHOD IS NOT SAFE MF20171121
//	public void clearDeviates() {
//		simulatedParameters.clear();
//		simulatedRandomEffects.clear();
//		simulatedResidualError.clear();
//		intervalLists.clear();
//		cruiseLineMap.clear();
//	}
	
	/**
	 * This method returns the residual error or the vector of residual errors associated with the subjectId.
	 * If the subject parameter is entered as null, the method assumes there is no need to store the simulated
	 * error terms in the simulatedResidualError map. This feature is useful if the residual error terms are 
	 * identically and independently distributed.
	 * @param subject a MonteCarloSimulationCompliantObject instance
	 * @param group an Enum that defines the group in case of different error term specifications
	 * @return a Matrix instance
	 */
	protected synchronized final Matrix getResidualErrorForThisSubject(MonteCarloSimulationCompliantObject subject, Enum<?> group) {
		if (isResidualVariabilityEnabled) {				// running in Monte Carlo mode
//			if (!rememberRandomDeviates) {
//				simulatedResidualError.clear();
//			}
			if (subject!= null && subject instanceof IndexableErrorTerm && defaultResidualError.get(group).getDistribution().isStructured()) {
				IndexableErrorTerm indexable = (IndexableErrorTerm) subject;
				GaussianErrorTermList list = getGaussianErrorTerms(subject);
				if (!list.getDistanceIndex().contains(indexable.getErrorTermIndex())) {
					list.add(new GaussianErrorTerm(indexable));
				}
				Matrix randomDeviate = defaultResidualError.get(group).getRandomDeviate(list);
				fireModelBasedSimulatorEvent(new REpiceaPredictorEvent(ModelBasedSimulatorEventProperty.RESIDUAL_ERROR_DEVIATE_JUST_GENERATED, null, new Object[]{subject, group, randomDeviate.getDeepClone()}, this));
				return randomDeviate; 
			} else {
				Matrix randomDeviate = defaultResidualError.get(group).getRandomDeviate();
				fireModelBasedSimulatorEvent(new REpiceaPredictorEvent(ModelBasedSimulatorEventProperty.RESIDUAL_ERROR_DEVIATE_JUST_GENERATED, null, new Object[]{subject, group, randomDeviate.getDeepClone()}, this));
				return randomDeviate;
			}
		} else {
			return defaultResidualError.get(group).getMean();
		}
	}
	
	protected synchronized final GaussianErrorTermList getGaussianErrorTerms(MonteCarloSimulationCompliantObject subject) {
		if (!doesThisSubjectHaveResidualErrorTerm(subject)) {		// the simulated parameters remain constant within the same Monte Carlo iteration
			simulatedResidualError.put(getSubjectPlusMonteCarloSpecificId(subject), new GaussianErrorTermList());
		}
		GaussianErrorTermList list =  simulatedResidualError.get(getSubjectPlusMonteCarloSpecificId(subject));
		return list;
	}
	
	protected final boolean doesThisSubjectHaveResidualErrorTerm(MonteCarloSimulationCompliantObject subject) {
		return simulatedResidualError.containsKey(getSubjectPlusMonteCarloSpecificId(subject));
	}
	
	
	/**
	 * This method returns the residual error under the assumption of iid and that the error is unique for all groups. 
	 * See the getResidualErrorForThisSubject method.
	 * @return a Matrix instance
	 */
	protected final Matrix getResidualError() {
		return getResidualErrorForThisSubject(null, ErrorTermGroup.Default);
	}

	protected void fireModelBasedSimulatorEvent(REpiceaPredictorEvent event) {
		for (REpiceaPredictorListener listener : listeners) {
			listener.modelBasedSimulatorDidThis(event);
		}
	}
	
//	protected void registerBlups(Matrix mean, Matrix variance, Matrix covariance, List<MonteCarloSimulationCompliantObject> subjectList) {
//		getParameterEstimates().registerBlups(mean, variance, covariance, subjectList);
//	}
	
	/**
	 * This method adds the listener instance to the list of listeners.
	 * @param listener a ModelBasedSimulatorListener listener
	 */
	public void addModelBasedSimulatorListener(REpiceaPredictorListener listener) {
		if (!listeners.contains(listener)) {
			listeners.add(listener);
		}
	}

	/**
	 * This method removes the listener instance from the list of listeners.
	 * @param listener a ModelBasedSimulatorListener listener
	 */
	public void removeModelBasedSimulatorListener(REpiceaPredictorListener listener) {
		listeners.remove(listener);
	}
	
//	/**
//	 * This method enables the recording of the random deviates. By default, this option is set to true.
//	 * It can be desirable to set this option to false when running large stochastic simulations.
//	 * @param rememberRandomDeviates a boolean
//	 */
//	public void setRememberRandomDeviates(boolean rememberRandomDeviates) {
//		this.rememberRandomDeviates = rememberRandomDeviates;
//	}

	protected boolean doBlupsExistForThisSubject(MonteCarloSimulationCompliantObject subject) {
		return getBlupsForThisSubject(subject) != null;
	}


	/**
	 * This method returns the blups for the subject or nothing if there is no blups for this subject
	 * @param subject a MonteCarloSimulationCompliantObject instance
	 * @return an Estimate instance or null
	 */
	protected Estimate<? extends StandardGaussianDistribution> getBlupsForThisSubject(MonteCarloSimulationCompliantObject subject) {
		if (blupsRandomEffects.containsKey(subject.getHierarchicalLevel().getName())) {
			if (blupsRandomEffects.get(subject.getHierarchicalLevel().getName()).containsKey(subject.getSubjectId())) {
				return blupsRandomEffects.get(subject.getHierarchicalLevel().getName()).get(subject.getSubjectId());
			}
		}
		return null;
	}

	protected final void setBlupsForThisSubject(MonteCarloSimulationCompliantObject subject, Estimate<? extends StandardGaussianDistribution> blups) {
		if (!blupsRandomEffects.containsKey(subject.getHierarchicalLevel().getName())) {
			blupsRandomEffects.put(subject.getHierarchicalLevel().getName(), new HashMap<String, Estimate<? extends StandardGaussianDistribution>>());
		}
		blupsRandomEffects.get(subject.getHierarchicalLevel().getName()).put(subject.getSubjectId(), blups);
		
		REpiceaPredictorEvent event = new REpiceaPredictorEvent(ModelBasedSimulatorEventProperty.BLUPS_JUST_SET, 
				null, 
				new Object[]{defaultRandomEffects.get(subject.getHierarchicalLevel().getName()), subject}, 
				this);
		fireModelBasedSimulatorEvent(event);
	}
	
	protected final void recordSubjectTestedForBlups(MonteCarloSimulationCompliantObject subject) {
		if (!subjectTestedForBlups.containsKey(subject.getHierarchicalLevel().getName())) {
			subjectTestedForBlups.put(subject.getHierarchicalLevel().getName(), new ArrayList<String>());
		}
		if (subjectTestedForBlups.get(subject.getHierarchicalLevel().getName()).contains(subject.getSubjectId())) {
			throw new InvalidParameterException("The subject has already been tested for blups!");
		} else {
			subjectTestedForBlups.get(subject.getHierarchicalLevel().getName()).add(subject.getSubjectId());
		}
	}

	protected final boolean hasSubjectBeenTestedForBlups(MonteCarloSimulationCompliantObject subject) {
		if (subjectTestedForBlups.containsKey(subject.getHierarchicalLevel().getName())) {
			return subjectTestedForBlups.get(subject.getHierarchicalLevel().getName()).contains(subject.getSubjectId());
		} else {
			return false;
		}
	}

	@Override
	public final boolean isStochastic() {
		return super.isStochastic() || isRandomEffectsVariabilityEnabled || isResidualVariabilityEnabled;
	}
	
}

