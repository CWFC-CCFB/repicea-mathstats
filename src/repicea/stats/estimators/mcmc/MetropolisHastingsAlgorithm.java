/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2021 Mathieu Fortin for Rouge Epicea.
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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;

import repicea.io.FormatField;
import repicea.io.javacsv.CSVField;
import repicea.io.javacsv.CSVWriter;
import repicea.math.Matrix;
import repicea.stats.StatisticalUtility;
import repicea.stats.data.DataSet;
import repicea.stats.distributions.GaussianDistribution;
import repicea.stats.estimates.MonteCarloEstimate;
import repicea.stats.estimators.AbstractEstimator;
import repicea.util.REpiceaLogManager;

/**
 * An implementation of the MCMC Metropolis-Hastings algorithm.
 * @author Mathieu Fortin - September 2021
 */
public class MetropolisHastingsAlgorithm  extends AbstractEstimator<MetropolisHastingsCompatibleModel> {
		
	private String loggerName;
	private String loggerPrefix;
	
	protected MetropolisHastingsParameters simParms;
	protected final MetropolisHastingsPriorHandler priors;
//	protected final MetropolisHastingsCompatibleModel model;
	private Matrix parameters;
	private Matrix parmsVarCov;
	protected double lpml;
	
	protected List<MetropolisHastingsSample> finalMetropolisHastingsSampleSelection;
	private boolean converged;
	protected int indexCorrelationParameter;
	private MonteCarloEstimate mcmcEstimate;

	public MetropolisHastingsAlgorithm(MetropolisHastingsCompatibleModel model, String loggerName, String loggerPrefix) {
		this(model);
		this.loggerName = loggerName;
		this.loggerPrefix = loggerPrefix;
	}
	
	public MetropolisHastingsAlgorithm(MetropolisHastingsCompatibleModel model) {
		super(model);
		simParms = new MetropolisHastingsParameters();
		priors = new MetropolisHastingsPriorHandler();
	}

	public void exportMetropolisHastingsSample(String filename) throws IOException {
		if (isConvergenceAchieved() && finalMetropolisHastingsSampleSelection != null) {
			CSVWriter writer = null;
			for (MetropolisHastingsSample sample : finalMetropolisHastingsSampleSelection) {
				if (writer == null) {
					writer = new CSVWriter(new File(filename), false);
					List<FormatField> fieldNames = new ArrayList<FormatField>();
					fieldNames.add(new CSVField("LLK"));
					for (int j = 1; j <= sample.parms.m_iRows; j++) {
						fieldNames.add(new CSVField("p" + j));
					}
					writer.setFields(fieldNames);
				}
				Object[] record = new Object[sample.parms.m_iRows + 1];
				record[0] = sample.llk;
				for (int j = 1; j <= sample.parms.m_iRows; j++) {
					record[j] = sample.parms.getValueAt(j - 1, 0);
				}
				writer.addRecord(record);
			}
			writer.close();
		}
	}
	
	public void setSimulationParameters(MetropolisHastingsParameters simParms) {
		if (simParms != null) {
			this.simParms = simParms;
		}
	}

	private String getLoggerName() {
		if (loggerName == null) {
			loggerName = getClass().getName();
		}
		return loggerName;
	}
	
	private MetropolisHastingsSample findFirstSetOfParameters(int desiredSize, boolean isForIntegral) {
		long startTime = System.currentTimeMillis();
		double llk = Double.NEGATIVE_INFINITY;
		List<MetropolisHastingsSample> myFirstList = new ArrayList<MetropolisHastingsSample>();
		while (myFirstList.size() < desiredSize) {
			Matrix parms = priors.getRandomRealization();
			llk = isForIntegral ? model.getLogLikelihood(parms) : model.getLogLikelihood(parms) + priors.getLogProbabilityDensityOfRandomEffects(parms) + priors.getLogProbabilityDensity(parms); // if isForIntegral then there is no need for the density of the parameters since the random realizations account for the distribution of the prior 
			if (llk > Double.NEGATIVE_INFINITY) {
				myFirstList.add(new MetropolisHastingsSample(parms, llk));
				if (myFirstList.size()%1000 == 0) {
					REpiceaLogManager.logMessage(getLoggerName(), Level.FINE, getLogMessagePrefix(), "Initial sample list has " + myFirstList.size() + " sets.");
				}
			}
		}
 		Collections.sort(myFirstList);
		MetropolisHastingsSample startingParms = myFirstList.get(myFirstList.size() - 1);
		REpiceaLogManager.logMessage(getLoggerName(), Level.FINE, getLogMessagePrefix(), "Time to find a first set of plausible parameters = " + (System.currentTimeMillis() - startTime) + " ms");
		REpiceaLogManager.logMessage(getLoggerName(), Level.FINE, getLogMessagePrefix(), "LLK = " + startingParms.llk + " - Parameters = " + startingParms.parms);
		return startingParms;
	}

	private String getLogMessagePrefix() {
		return loggerPrefix == null ? "" : loggerPrefix;
	}

	public MetropolisHastingsParameters getSimulationParameters() {
		return simParms;
	}

//	public double getMarginalLogLikelihood() {
//		return lnProbY;
//	}
//
//	public double getMarginalLogLikelihood2() {
//		return lnProbY2;
//	}

	/**
	 * Return the final parameter estimates. <br>
	 * <br>
	 * Convergence must be achieved. If the parameters member has not been set, it is then
	 * set on the fly to avoid recalculating the mean from the MonteCarloEstimate every time.
	 * @return a Matrix
	 */
	public Matrix getFinalParameterEstimates() {
		if (isConvergenceAchieved()) {
			if (parameters == null) 
				parameters = mcmcEstimate.getMean();
			return parameters;
		} else {
			return null;
		}
	}
	
	/**
	 * Return the estimated variance-covariance of the final parameter estimates. <br>
	 * <br>
	 * Convergence must be achieved. If the variance-covariance member has not been set, it is then
	 * set on the fly to avoid recalculating it from the MonteCarloEstimate every time.
	 * @return a Matrix
	 */
	public Matrix getParameterCovarianceMatrix() {
		if (isConvergenceAchieved()) {
			if (parmsVarCov == null) 
				parmsVarCov = mcmcEstimate.getVariance();
			return parmsVarCov;
		} else {
			return null;
		}
	}
	
	@Override
	public boolean isConvergenceAchieved() {
		return converged;
	}
	
	public MetropolisHastingsPriorHandler getPriorHandler() {
		return priors;
	}
	
	/**
	 * Implement the Metropolis-Hastings algorithm.
	 * @param nbRealizations number of samples in the chain before removing burn in and selected one sample every x samples
	 * @param nbBurnIn number of samples to discard at the beginning of the chain
	 * @param nbInternalIter maximum number of realizations to find the next acceptable sample of the chain
	 * @param metropolisHastingsSample A list of MetaModelMetropolisHastingsSample instance that represents the chain
	 * @param gaussDist the sampling distribution
	 * @return a boolean
	 */
	private boolean generateMetropolisSample(List<MetropolisHastingsSample> metropolisHastingsSample, GaussianDistribution gaussDist) {
		long startTime = System.currentTimeMillis();
		Matrix newParms = null;
		double llk = 0d;
		boolean completed = true;
		int trials = 0;
		int successes = 0;
		double acceptanceRatio; 
		for (int i = 0; i < simParms.nbRealizations - 1; i++) { // Metropolis-Hasting  -1 : the starting parameters are considered as the first realization
			gaussDist.setMean(metropolisHastingsSample.get(metropolisHastingsSample.size() - 1).parms);
			if (i > 0 && i < simParms.nbBurnIn && i%1000 == 0) {
				acceptanceRatio = ((double) successes) / trials;
				REpiceaLogManager.logMessage(getLoggerName(), Level.FINE, getLogMessagePrefix(), "After " + i + " realizations, the acceptance rate is " + acceptanceRatio);
				if (acceptanceRatio > 0.40) {	// we aim at having an acceptance rate slightly larger than 0.3 because it will decrease as the chain reaches its steady state
					gaussDist.setVariance(gaussDist.getVariance().scalarMultiply(1.2*1.2));
				} else if (acceptanceRatio < 0.30) {
					gaussDist.setVariance(gaussDist.getVariance().scalarMultiply(0.8*0.8));
				}
				successes = 0;
				trials = 0;
			}
			if (i%10000 == 0 && i > simParms.nbBurnIn) {
				acceptanceRatio = ((double) successes) / trials;
				REpiceaLogManager.logMessage(getLoggerName(), Level.FINE, getLogMessagePrefix(), "Processing realization " + i + " / " + simParms.nbRealizations + "; " + acceptanceRatio);
			}
			boolean accepted = false;
			int innerIter = 0;
			
			while (!accepted && innerIter < simParms.nbInternalIter) {
				newParms = gaussDist.getRandomRealization();
				double parmsPriorLogDensity = priors.getLogProbabilityDensity(newParms);
				if (parmsPriorLogDensity > Double.NEGATIVE_INFINITY) {
					llk = model.getLogLikelihood(newParms) + 
							priors.getLogProbabilityDensityOfRandomEffects(newParms) + 
							parmsPriorLogDensity;
					double ratio = Math.exp(llk - metropolisHastingsSample.get(metropolisHastingsSample.size() - 1).llk);
					accepted = StatisticalUtility.getRandom().nextDouble() < ratio;
					trials++;
					if (accepted) {
						successes++;
					}
				}
				innerIter++;
			}
			if (innerIter >= simParms.nbInternalIter && !accepted) {
				REpiceaLogManager.logMessage(getLoggerName(), Level.SEVERE,  getLogMessagePrefix(), "Stopping after " + i + " realization");
				completed = false;
				break;
			} else {
				metropolisHastingsSample.add(new MetropolisHastingsSample(newParms, llk));  // new set of parameters is recorded
				if (metropolisHastingsSample.size()%100 == 0) {
					REpiceaLogManager.logMessage(getLoggerName(), Level.FINEST, getLogMessagePrefix(), metropolisHastingsSample.get(metropolisHastingsSample.size() - 1));
				}
			}
		}
		
		if (completed) {
			acceptanceRatio = ((double) successes) / trials;
			REpiceaLogManager.logMessage(getLoggerName(), Level.INFO, getLogMessagePrefix(), "Time to obtain " + metropolisHastingsSample.size() + " samples = " + (System.currentTimeMillis() - startTime) + " ms");
			REpiceaLogManager.logMessage(getLoggerName(), Level.INFO, getLogMessagePrefix(), "Acceptance ratio = " + acceptanceRatio);
		} 
		return completed;
	}

	private void resetSuccessAndTrialMaps(GaussianDistribution dist, 
			Map<Integer, Integer> trials, 
			Map<Integer, Integer> successes) {
		trials.clear();
		successes.clear();
		for (int i = 0; i < dist.getMean().m_iRows; i++) {
			trials.put(i, 0);
			successes.put(i, 0);
		}
	}
	
	private Matrix computeSuccessRates(Map<Integer, Integer> trials, Map<Integer, Integer> successes) {
		Matrix ratios = new Matrix(trials.size(), 1);
		for (int i = 0; i < trials.size(); i++) {
			double ratio = ((double) successes.get(i)) / trials.get(i);
			ratios.setValueAt(i, 0, ratio);
		}
		return ratios;
	}
	
	/**
	 * Implement Gibbs sampling in a preliminary stage to balance the variance of the sampler.
	 * @param firstSample the MetaModelMetropolisHastingsSample instance that was found through random sampling
	 * @param sampler the sampling distribution
	 * @return a boolean
	 */
	private boolean balanceVariance(MetropolisHastingsSample firstSample, GaussianDistribution sampler) {
		long startTime = System.currentTimeMillis();
		List<MetropolisHastingsSample> initSample = new ArrayList<MetropolisHastingsSample>();
		initSample.add(firstSample);
		Matrix newParms = null;
		double llk = 0d;
		boolean completed = true;
		Matrix acceptanceRatios; 
		Map<Integer, Integer> trialMap = new HashMap<Integer, Integer>();
		Map<Integer, Integer> successMap = new HashMap<Integer, Integer>();
		resetSuccessAndTrialMaps(sampler, trialMap, successMap);
		double targetAcceptance = 0.5; // MF2021-11-01 This number does not matter much in absolute value. It just makes sure that the acceptance rate is balanced across the parameters.
		for (int i = 0; i < simParms.nbBurnIn - 1; i++) { // Metropolis-Hasting  -1 : the starting parameters are considered as the first realization
			Matrix originalParms = initSample.get(initSample.size() - 1).parms.getDeepClone();
			sampler.setMean(originalParms);
			if (i > 0 && i < simParms.nbBurnIn && i%1000 == 0) {
				acceptanceRatios = this.computeSuccessRates(trialMap, successMap);
				REpiceaLogManager.logMessage(getLoggerName(), Level.FINE, getLogMessagePrefix(), "After " + i + " realizations, the acceptance rates are " + acceptanceRatios);
				for (int j = 0; j < acceptanceRatios.m_iRows; j++) {
					double currentRatio = acceptanceRatios.getValueAt(j, 0);
					if (currentRatio > targetAcceptance + .05) {	// then we must increase the CoefVar
						Matrix variance = sampler.getVariance();
						variance.setValueAt(j, j, variance.getValueAt(j, j) * 1.2 * 1.2);
					} else if (currentRatio < targetAcceptance - .05) {
						Matrix variance = sampler.getVariance();
						variance.setValueAt(j, j, variance.getValueAt(j, j) * 0.8 * 0.8);
					}
				}
				resetSuccessAndTrialMaps(sampler, trialMap, successMap);
			}
			boolean accepted = false;
			int innerIter = 0;

			int j = 0;
			while (j < originalParms.m_iRows && innerIter < simParms.nbInternalIter) {
				double originalValue = originalParms.getValueAt(j, 0);
				double newValue = getNewParms(sampler, j);
				originalParms.setValueAt(j, 0, newValue);
				double parmsPriorLogDensity = priors.getLogProbabilityDensity(originalParms);
				if (parmsPriorLogDensity > Double.NEGATIVE_INFINITY) {
					llk = model.getLogLikelihood(originalParms) + 
							priors.getLogProbabilityDensityOfRandomEffects(originalParms) +
							parmsPriorLogDensity;
					double ratio = Math.exp(llk - initSample.get(initSample.size() - 1).llk);
					accepted = StatisticalUtility.getRandom().nextDouble() < ratio;
					trialMap.put(j, trialMap.get(j) + 1);
					if (accepted) {
						successMap.put(j, successMap.get(j) + 1);
						j++;
						accepted = false;
						innerIter = 0;
					} else {
						originalParms.setValueAt(j, 0, originalValue);	// we put the old value back into the vector of parameters
					}
				} else {
					originalParms.setValueAt(j, 0, originalValue);	// we put the old value back into the vector of parameters
				}
				innerIter++;
			}
			newParms = originalParms;
			if (innerIter >= simParms.nbInternalIter && !accepted) {
				REpiceaLogManager.logMessage(getLoggerName(), Level.SEVERE,  getLogMessagePrefix(), "Stopping after " + i + " realization");
				completed = false;
				break;
			} else {
				initSample.add(new MetropolisHastingsSample(newParms, llk));  // new set of parameters is recorded
				if (initSample.size()%100 == 0) {
					REpiceaLogManager.logMessage(getLoggerName(), Level.FINEST, getLogMessagePrefix(), initSample.get(initSample.size() - 1));
				}
			}
		}
		
		if (completed) {
			acceptanceRatios = computeSuccessRates(trialMap, successMap);
			REpiceaLogManager.logMessage(getLoggerName(), Level.FINE, getLogMessagePrefix(), "Time to balance the variance of the sample: " + (System.currentTimeMillis() - startTime) + " ms");
			REpiceaLogManager.logMessage(getLoggerName(), Level.FINE, getLogMessagePrefix(), "Acceptance ratio = " + acceptanceRatios);
		} 
		return completed;
	}

	private double getNewParms(GaussianDistribution dist, int i) {
		double variance = dist.getVariance().getValueAt(i, i);
		double mean = dist.getMean().getValueAt(i, 0);
		double newValue = mean + StatisticalUtility.getRandom().nextGaussian() * Math.sqrt(variance);
		return newValue;
	}
	
	private List<MetropolisHastingsSample> retrieveFinalSample(List<MetropolisHastingsSample> metropolisHastingsSample) {
		List<MetropolisHastingsSample> finalMetropolisHastingsGibbsSample = new ArrayList<MetropolisHastingsSample>();
		REpiceaLogManager.logMessage(getLoggerName(), Level.FINE, getLogMessagePrefix(), "Discarding " + simParms.nbBurnIn + " samples as burn in.");
		for (int i = simParms.nbBurnIn; i < metropolisHastingsSample.size(); i+= simParms.oneEach) {
			finalMetropolisHastingsGibbsSample.add(metropolisHastingsSample.get(i));
		}
		REpiceaLogManager.logMessage(getLoggerName(), Level.FINE, getLogMessagePrefix(), "Selecting one every " + simParms.oneEach + " samples as final selection.");
		return finalMetropolisHastingsGibbsSample;
	}

	private void reset() {
		parameters = null;
		parmsVarCov = null;
		mcmcEstimate = null;
//		lnProbY = 0;
		finalMetropolisHastingsSampleSelection = null;
		converged = false;
	}
	
	@Override
	public boolean doEstimation() {
		reset();
		double coefVar = 0.01;
		try {
			GaussianDistribution samplingDist = model.getStartingParmEst(coefVar);
			List<MetropolisHastingsSample> mhSample = new ArrayList<MetropolisHastingsSample>();
			MetropolisHastingsSample firstSet = findFirstSetOfParameters(simParms.nbInitialGrid, false);	// false: not for integration
			mhSample.add(firstSet); // first valid sample
			boolean completed = balanceVariance(firstSet, samplingDist);
			if (completed) {
				completed = generateMetropolisSample(mhSample, samplingDist);
				if (completed) {
					finalMetropolisHastingsSampleSelection = retrieveFinalSample(mhSample);
					mcmcEstimate = new MonteCarloEstimate();
					for (MetropolisHastingsSample sample : finalMetropolisHastingsSampleSelection) {
						mcmcEstimate.addRealization(sample.parms);
					}

					List<MetropolisHastingsSample> tempSample = new ArrayList<MetropolisHastingsSample>();
					tempSample.addAll(finalMetropolisHastingsSampleSelection);
					Collections.sort(tempSample);
					this.lpml = calculateLogPseudomarginalLikelihood();
					REpiceaLogManager.logMessage(getLoggerName(), Level.FINE, getLogMessagePrefix(), "Final sample had " + finalMetropolisHastingsSampleSelection.size() + " sets of parameters.");
					converged = true;
				}
			}
		} catch (Exception e1) {
			e1.printStackTrace();
			converged = false;
		} 
		return converged;
	}
	

	private double calculateLogPseudomarginalLikelihood() {
		int nbSubjects = model.getNbSubjects();
		double lpml = 0;
		for (int i = 0; i < nbSubjects; i++) {
			double sum = 0;
			for (MetropolisHastingsSample s : finalMetropolisHastingsSampleSelection) {
				double lk_i = model.getLikelihoodOfThisSubject(s.parms, i);
				sum += 1d / lk_i;
			}
			sum /= finalMetropolisHastingsSampleSelection.size();
			double cpo_i = 1d / sum;
			lpml += Math.log(cpo_i);
		}
//		lpml /= nbSubjects;			// This division is a typo in de la Cruz et al. (2016). Validated with Hoque (2017) and Ibrahim et al. (2001 p.228).
		return lpml;
	}
	
	public double getLogPseudomarginalLikelihood() {
		if (isConvergenceAchieved()) 
			return lpml;
		else return Double.NaN;
	}

	@Override
	public MonteCarloEstimate getParameterEstimates() {
		return isConvergenceAchieved() ? mcmcEstimate : null;
	}
	
	@Override
	public DataSet getConvergenceStatusReport() {
		DataSet ds = super.getConvergenceStatusReport();
		Object[] record = new Object[2];
		record[0] = "Log Pseudomarginal Likelihood";
		record[1] = lpml;
		ds.addObservation(record);
		return ds;
	}


	
	
}
