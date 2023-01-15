/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2020 Mathieu Fortin for Rouge-Epicea
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

import java.lang.reflect.Constructor;
import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.List;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.distributions.EmpiricalDistribution;
import repicea.stats.distributions.UnknownDistribution;
import repicea.stats.distributions.utility.GaussianUtility;
import repicea.stats.sampling.PopulationUnitWithEqualInclusionProbability;
import repicea.stats.sampling.PopulationUnitWithUnequalInclusionProbability;

/**
 * This class implements the bootstrap estimator of the total as developed by Fortin et al. 
 * (2018) in the context of hybrid inference. More specifically, it applies when (i) the 
 * variable of interest has not been observed but predicted by a model and (ii) the covariates 
 * are not censused but only observed in a sample of the population. The estimator takes into 
 * account the variance that stems from the model as well as that of the sampling design. 
 * IMPORTANT: The model must  benefit from a full stochastic implementation.
 * @author Mathieu Fortin - May 2018
 * @see <a href=https://academic.oup.com/forestry/article/91/3/354/4647707>
 * Fortin, M., Manso, R., and Schneider, R. 2018. Parametric bootstrap estimators for hybrid 
 * inference in forest inventories. Forestry 91(3): 354-365. </a>
 */
@SuppressWarnings("serial")
public final class BootstrapHybridPointEstimate extends Estimate<UnknownDistribution> implements NumberOfRealizationsProvider {

	private static boolean IsVarianceCorrectionEnabled = true; 	// default value
	
	public static enum VarianceEstimatorImplementation {
		/**
		 * The unbiased estimator in Fortin et al. (2017).
		 */
		Corrected,
		/**
		 * A less biased estimator than the multiple imputation. 
		 */
		LessBiased,
		/**
		 * The estimator in Rubin (1987).
		 */
		RegularMultipleImputation;
	}

	/**
	 * An inner class that contains the corrected variance of the BootstrapHybridPointEstimate instance.
	 * @author Mathieu Fortin - May 2018
	 */
	class VarianceEstimate extends SimpleEstimate {
		
		final int numberOfRealizations;
		final SymmetricMatrix varMean;
		final SymmetricMatrix meanVar;
		final SymmetricMatrix designVarianceOfMeanRealizedY;
		final SymmetricMatrix varianceBiasCorrection;
		final boolean isCalculable;
		
		/**
		 * Private constructor for the different variance estimator implementation. <br>
		 * <br>
		 * There are three implementations: Corrected, LessBiased, RegularMultipleImputation and None. If the argument
		 * designVarianceOfMeanRealizedY is non null, the constructor tries to corrected implementation. If this 
		 * implementation fails, that is if it provides negative variance estimates, then the constructor relies
		 * on the less biased implementation. If designVarianceOfMeanRealizedY is set to null then the implementation
		 * is the regular multiple imputation which is the most biased. <br>
		 * <br>
		 * If any of the varMean or meanVar argument is null, then the implementation is set to None.
		 * 
		 * @param pointEstimate the point estimate 
		 * @param varMean the variance of the bootstrapped point estimates 
		 * @param meanVar the mean of the bootstrapped variances
		 * @param designVarianceOfMeanRealizedY the design variance based on the mean realizations of y. 
		 * @param rowIndex a list of strings that represent the row indices
		 */
		VarianceEstimate(int numberOfRealizations,
				Matrix pointEstimate, 
				SymmetricMatrix varMean, 
				SymmetricMatrix meanVar, 
				SymmetricMatrix designVarianceOfMeanRealizedY, 
				List<String> rowIndex) {
			super(pointEstimate, null); 
			if (pointEstimate == null) {
				throw new InvalidParameterException("The pointEstimate argument must be non null!");
			}
			this.numberOfRealizations = numberOfRealizations;
			this.varMean = varMean != null ? varMean.getDeepClone() : null;
			this.meanVar = varMean != null ? meanVar.getDeepClone() : null;
			this.designVarianceOfMeanRealizedY = designVarianceOfMeanRealizedY != null ? designVarianceOfMeanRealizedY.getDeepClone() : null;
			if (this.varMean!= null && this.meanVar != null && this.designVarianceOfMeanRealizedY != null) {
				isCalculable = true;
			} else {
				isCalculable = false;
			}
			if (this.varMean != null && this.meanVar != null && this.designVarianceOfMeanRealizedY != null) {
				SymmetricMatrix varMeanPlusDesignVarianceOfMeanRealizedY = (SymmetricMatrix) this.varMean.scalarMultiply((this.numberOfRealizations - 1d) / this.numberOfRealizations).add(this.designVarianceOfMeanRealizedY);		// factor (n-1)/n in order to get the sum of square divided by n and not by n-1
				SymmetricMatrix modelRelatedVariance = (SymmetricMatrix) varMeanPlusDesignVarianceOfMeanRealizedY.subtract(this.meanVar);
				
				if (modelRelatedVariance.diagonalVector().anyElementSmallerOrEqualTo(0d) || !BootstrapHybridPointEstimate.IsVarianceCorrectionEnabled) { // means that the corrected variance estimator is inconsistent or it has been overriden
					this.varianceBiasCorrection = null;
				} else {
					Matrix denominator = pointEstimate.multiply(pointEstimate.transpose()).subtract(this.designVarianceOfMeanRealizedY);
					Matrix numerator = modelRelatedVariance.elementWiseMultiply(this.designVarianceOfMeanRealizedY);
					varianceBiasCorrection = SymmetricMatrix.convertToSymmetricIfPossible(numerator.elementWiseDivide(denominator).scalarMultiply(-1d));
				}
			} else {
				this.varianceBiasCorrection = null;
			}
			setRowIndex(rowIndex);  
		}

		@Override
		protected SymmetricMatrix getVarianceFromDistribution() {
			if (isCalculable) {
				Matrix varMeanPlusDesignVarianceOfMeanRealizedY = varMean.
						scalarMultiply((numberOfRealizations - 1d) / numberOfRealizations).
						add(designVarianceOfMeanRealizedY);		// factor (n-1)/n in order to get the sum of square divided by n and not by n-1
				Matrix modelRelatedVariance = (SymmetricMatrix) varMeanPlusDesignVarianceOfMeanRealizedY.subtract(meanVar);
				switch(BootstrapHybridPointEstimate.this.vei) {
				case Corrected:
					return (SymmetricMatrix) modelRelatedVariance.add(designVarianceOfMeanRealizedY);
				case LessBiased:
					return (SymmetricMatrix) varMeanPlusDesignVarianceOfMeanRealizedY;
				case RegularMultipleImputation:
					return (SymmetricMatrix) varMean.scalarMultiply((numberOfRealizations + 1d) / numberOfRealizations).add(meanVar);	// the (n + 1)/n factor comes from Rubin 1987 p.76;
				default:
					throw new InvalidParameterException("This variance estimator implementation is unknown: " + BootstrapHybridPointEstimate.this.vei.name());
				}
			} else {
				return null;
			}
		}
		
		protected Matrix getVarianceBiasCorrection() {return varianceBiasCorrection;}

		
		/**
		 * Provide the model-related variance, which includes the variance bias correction if the implementation
		 * is the corrected one. 
		 * @return a Matrix instance
		 */
		Matrix getModelRelatedVariance() {
			if (isCalculable) {
				switch(BootstrapHybridPointEstimate.this.vei) {
				case Corrected:
					return varMean.scalarMultiply((numberOfRealizations - 1d) / numberOfRealizations).add(designVarianceOfMeanRealizedY).subtract(meanVar);
				case LessBiased:
					return varMean.scalarMultiply((numberOfRealizations - 1d) / numberOfRealizations);
				case RegularMultipleImputation:
					return varMean.scalarMultiply((numberOfRealizations + 1d) / numberOfRealizations);
				default:
					throw new InvalidParameterException("This variance estimator implementation is unknown: " + BootstrapHybridPointEstimate.this.vei.name());
				}
			} else {
				return null;
			}
		}

		
		/**
		 * Provide the model-related variance, without the variance bias correction. <br>
		 * <br>
		 * This method assumes that the estimator implementation is the corrected one. Otherwise it returns null.
		 * 
		 * @return a Matrix instance or null if the estimator implementation is not EstimatorImplementation.Corrected.
		 */
		public Matrix getNetModelRelatedVariance() {
			if (getVarianceBiasCorrection() != null) {
				return getModelRelatedVariance().subtract(getVarianceBiasCorrection());
			} else {
				return null;
			}
		}
		
		/**
		 * This method returns the estimate of the sampling-related variance.
		 * @return a Matrix instance
		 */
		Matrix getSamplingRelatedVariance() {
			if (isCalculable) {
				switch(BootstrapHybridPointEstimate.this.vei) {
				case Corrected:
					return designVarianceOfMeanRealizedY;
				case LessBiased:
					return designVarianceOfMeanRealizedY;
				case RegularMultipleImputation:
					return meanVar;
				default:
					throw new InvalidParameterException("This variance estimator implementation is unknown: " + BootstrapHybridPointEstimate.this.vei.name());
				}
			} else {
				return null;
			}
		}
		
	}
	
	
	private final List<PointEstimate<?>> estimates;
	private VarianceEstimate varianceEstimate;
	private VarianceEstimatorImplementation vei;

	/**
	 * Constructor.
	 */
	public BootstrapHybridPointEstimate() {
		super(new UnknownDistribution());
		estimates = new ArrayList<PointEstimate<?>>();
		vei = VarianceEstimatorImplementation.Corrected; // default value
	}

//	/**
//	 * Allow to enable or disable the variance correction. If the variance correction is disabled,
//	 * then the less biased variance estimator is preferred. By default, the variance correction is enabled.
//	 * @param bool a boolean
//	 */
//	public static void setVarianceCorrectionEnabled(boolean bool) {	
//		BootstrapHybridPointEstimate.IsVarianceCorrectionEnabled = bool;
//	}

	/**
	 * Set the variance estimator implementation. Can be corrected, less biased or regular multiple imputation (or eventually none). 
	 * @param vei a VarianceEstimatorImplementation enum
	 */
	public void setVarianceEstimatorImplementation(VarianceEstimatorImplementation vei) {
		if (vei != null) {
			this.vei = vei;
		}
	}

	/**
	 * Return the current variance implementation. 
	 * @return a VarianceEstimatorImplementation enum
	 */
	public VarianceEstimatorImplementation getVarianceEstimatorImplementation() {
		return vei;
	}

	
	
	
	/**
	 * Return true if the variance correction is enabled (default) or false otherwise.
	 * @return a boolean
	 */
	public static boolean isVarianceCorrectionEnabled() {	
		return BootstrapHybridPointEstimate.IsVarianceCorrectionEnabled;
	}
	
	/**
	 * This method adds a realization of the point estimate. The compatibility of 
	 * the instance with previously added instances is checked. If the check fails
	 * an InvalidParameterException is thrown.
	 * @param estimate a PointEstimate instance
	 */
	public void addPointEstimate(PointEstimate<?> estimate) {
		if (estimates.isEmpty() || estimates.get(0).isMergeableEstimate(estimate)) {
			estimates.add(estimate);
			varianceEstimate = null;	// make sure to reset the variance estimate
		} else {
			throw new InvalidParameterException("The point estimate is not compatible with the previous estimates!");
		}
	}
	
	/**
	 * Provide the model-related variance, which includes the variance bias correction if the implementation
	 * is the corrected one. 
	 * @return a Matrix instance
	 */
	public Matrix getModelRelatedVariance() {
		return getVarianceEstimate().getModelRelatedVariance();
	}
	
	/**
	 * This method returns the estimate of the sampling-related variance.
	 * @return a Matrix instance
	 */
	public Matrix getSamplingRelatedVariance() {
		return getVarianceEstimate().getSamplingRelatedVariance();
	}
	
	/*
	 * For test only.
	 */
	protected Matrix getVarianceBiasCorrection() {
		return getVarianceEstimate().getVarianceBiasCorrection();
	}
	
	
	@Override
	protected boolean isMergeableEstimate(Estimate<?> estimate) {
		if (estimate instanceof BootstrapHybridPointEstimate) {
			if (((BootstrapHybridPointEstimate) estimate).getNumberOfRealizations() == getNumberOfRealizations()) {
				if (estimates.get(0).isMergeableEstimate(((BootstrapHybridPointEstimate) estimate).estimates.get(0)));
				return true;
			};
		}
		return false;
	}
	
	/**
	 * This method is useful for large simulations. The simulation can be run in different threads or batches
	 * and the resulting BootstrapHybridPointEstimate can then be combined through this method.
	 * @param estimate a BootstrapHybridPointEstimate instance that relies on the same PointEstimate class in the estimates member
	 */
	public void appendBootstrapHybridEstimate(BootstrapHybridPointEstimate estimate) {
		for (PointEstimate<?> pointEstimate : estimate.estimates) {
			addPointEstimate(pointEstimate);
		}
	}

	@Override
	protected Matrix getMeanFromDistribution() {
		Matrix mean = null;
		for (PointEstimate<?> estimate : estimates) {
			if (mean == null) {
				mean = estimate.getMean();
			} else {
				mean = mean.add(estimate.getMean());
			}
		}
		mean = mean.scalarMultiply(1d / estimates.size());
		return mean;
	}
	
	
	

//	/**
//	 * This method returns the uncorrected variance of the total estimate. 
//	 * This estimator is based on the law of total variance. It tends to overestimate 
//	 * the true variance. This method is deprecated and the getCorrectedVariance method
//	 * should be used in place.
//	 * @return a Matrix
//	 * @see <a href=https://academic.oup.com/forestry/article/91/3/354/4647707>
//	 * Fortin, M., Manso, R., and Schneider, R. 2018. Parametric bootstrap estimators for hybrid 
//	 * inference in forest inventories. Forestry 91(3): 354-365. </a>
//	 */
//	@Deprecated
//	public final VarianceEstimate getUncorrectedVariance() {
//		MonteCarloEstimate variance = new MonteCarloEstimate();
//		MonteCarloEstimate mean = new MonteCarloEstimate();
//		for (PointEstimate<?> estimate : estimates) {
//			mean.addRealization(estimate.getMean());
//			variance.addRealization(estimate.getVariance());
//		}
//		return new VariancePointEstimate(mean.getNumberOfRealizations(),
//				mean.getMean(), 
//				mean.getVariance(), 
//				variance.getMean(),
//				null,
//				rowIndex);
//	}

	
	/**
	 * This method calculates the corrected variance of the total estimate. 
	 * This estimator is theoretically unbiased. 
	 * @return a VariancePointEstimate
	 * 
	 * @see <a href=https://academic.oup.com/forestry/article/91/3/354/4647707>
	 * Fortin, M., Manso, R., and Schneider, R. 2018. Parametric bootstrap estimators for hybrid 
	 * inference in forest inventories. Forestry 91(3): 354-365. </a>
	 */
	private VarianceEstimate calculateVarianceEstimate() {
		if (getNumberOfRealizations() > 1) {
			MonteCarloEstimate variance = new MonteCarloEstimate();
			MonteCarloEstimate mean = new MonteCarloEstimate();
			int sampleSize = estimates.get(0).getObservations().size();
			EmpiricalDistribution observationMeans = new EmpiricalDistribution();
			int nbElementsPerObs = 0;
			for (PointEstimate<?> estimate : estimates) {
				if (nbElementsPerObs == 0) {
					nbElementsPerObs = estimate.getNumberOfElementsPerObservation();
				}
				mean.addRealization(estimate.getMean());
				variance.addRealization(estimate.getVariance());
				observationMeans.addRealization(estimate.getObservationMatrix());
			}
			
			PointEstimate<?> meanEstimate; 
			try {
				if (estimates.get(0).isPopulationSizeKnown()) {
					double populationSize = estimates.get(0).getPopulationSize();
					Constructor<?> cons = estimates.get(0).getClass().getConstructor(double.class);
					meanEstimate = (PointEstimate<?>) cons.newInstance(populationSize);
				} else {
					meanEstimate = estimates.get(0).getClass().newInstance();
				}
				
				Matrix observationMeanMatrix = observationMeans.getMean();
				List<String> sampleIds = estimates.get(0).getSampleIds();
				for (int i = 0; i < sampleSize; i++) {
					String sampleId = sampleIds.get(i);
					Matrix meanForThisI = observationMeanMatrix.getSubMatrix(i, i, 0, nbElementsPerObs - 1).transpose();		// need to transpose because it is a row vector
					if (meanEstimate instanceof PopulationTotalEstimate) {
						PopulationUnitWithUnequalInclusionProbability popUnit = new PopulationUnitWithUnequalInclusionProbability(sampleId, meanForThisI, ((PopulationTotalEstimate) estimates.get(0)).getObservations().get(sampleId).getInclusionProbability());
						((PopulationTotalEstimate) meanEstimate).addObservation(popUnit);
					} else {
						PopulationUnitWithEqualInclusionProbability popUnit = new PopulationUnitWithEqualInclusionProbability(sampleId, meanForThisI);
						((PopulationMeanEstimate) meanEstimate).addObservation(popUnit);
					}
				}
						
				VarianceEstimate varEst = new VarianceEstimate(getNumberOfRealizations(),
						getMean(),
						mean.getVariance(),
					    SymmetricMatrix.convertToSymmetricIfPossible(variance.getMean()), 
						meanEstimate.getVariance(), 
						rowIndex);
				return varEst;
			} catch (Exception e) {
				throw new InvalidParameterException("An error occured while instantiating the correct PointEstimate class!");
			}
		} else {
			System.out.println("The variance of the hybrid point estimate cannot be calculated because there is not enough realizations!");
			return new VarianceEstimate(getNumberOfRealizations(), getMean(), null, null, null, rowIndex);
		}
	}
	
	private VarianceEstimate getVarianceEstimate() {
		if (varianceEstimate == null) {
			varianceEstimate = calculateVarianceEstimate();
		}
		return varianceEstimate;
	}
	

	@Override
	protected SymmetricMatrix getVarianceFromDistribution() {
		return getVarianceEstimate().getVariance();
	}

	@Override
	public int getNumberOfRealizations() {
		return estimates.size();
	}
	
	protected Matrix getQuantileForProbability(double probability) {
		Matrix stdDev = getVariance().diagonalVector().elementWisePower(.5); 
		double quantile = GaussianUtility.getQuantile(probability);
		return getMean().add(stdDev.scalarMultiply(quantile));
	}

	@Override
	public ConfidenceInterval getConfidenceIntervalBounds(double oneMinusAlpha) {
		Matrix lowerBoundValue = getQuantileForProbability(.5 * (1d - oneMinusAlpha));
		Matrix upperBoundValue = getQuantileForProbability(1d - .5 * (1d - oneMinusAlpha));
		return new ConfidenceInterval(lowerBoundValue, upperBoundValue, oneMinusAlpha);
	}

	/**
	 * Returns the sample size the individual point estimates are based on or -1 if 
	 * no point estimate is available at the moment.
	 * @return an integer
	 */
	public int getSampleSize() {
		if (estimates == null || estimates.isEmpty()) {
			return -1;
		} else {
			return estimates.get(0).getObservations().size();
		}
	}
	
	
	
	
	/**
	 * This method returns a MonteCarloEstimate instance that results from the subtraction of two 
	 * MonteCarloEstimate instances with the same number of realizations. 
	 * @param estimate2 the estimate that is subtracted from this estimate
	 * @return a BootstrapHybridPointEstimate instance
	 */
	protected BootstrapHybridPointEstimate subtract(BootstrapHybridPointEstimate estimate2) {
		if (getNumberOfRealizations() != estimate2.getNumberOfRealizations()) {
			throw new InvalidParameterException("The number of realizations is not consistent!");
		}
		BootstrapHybridPointEstimate outputEstimate = new BootstrapHybridPointEstimate();
		for (int i = 0; i < getNumberOfRealizations(); i++) {
			outputEstimate.addPointEstimate(estimates.get(i).subtract(estimate2.estimates.get(i)));
		}
		return outputEstimate;
	}
	
	/**
	 * This method returns a MonteCarloEstimate instance that results from the sum of two 
	 * MonteCarloEstimate instances with the same number of realizations. 
	 * @param estimate2 the estimate that is added to this estimate
	 * @return a BootstrapHybridPointEstimate instance
	 */
	protected BootstrapHybridPointEstimate add(BootstrapHybridPointEstimate estimate2) {
		if (getNumberOfRealizations() != estimate2.getNumberOfRealizations()) {
			throw new InvalidParameterException("The number of realizations is not consistent!");
		}
		BootstrapHybridPointEstimate outputEstimate = new BootstrapHybridPointEstimate();
		for (int i = 0; i < getNumberOfRealizations(); i++) {
			outputEstimate.addPointEstimate(estimates.get(i).add(estimate2.estimates.get(i)));
		}
		return outputEstimate;
	}

	/**
	 * This method returns a MonteCarloEstimate instance that results from the product of original 
	 * MonteCarloEstimate instance and a scalar. 
	 * @param scalar the multiplication factor
	 * @return a BootstrapHybridPointEstimate instance
	 */
	protected BootstrapHybridPointEstimate multiply(double scalar) {
		BootstrapHybridPointEstimate outputEstimate = new BootstrapHybridPointEstimate();
		for (int i = 0; i < getNumberOfRealizations(); i++) {
			outputEstimate.addPointEstimate(estimates.get(i).multiply(scalar));
		}
		return outputEstimate;
	}

	@Override
	public Estimate<?> getDifferenceEstimate(Estimate<?> estimate2) {
		if (this.isMergeableEstimate(estimate2)) {
			return subtract((BootstrapHybridPointEstimate) estimate2);
		} else {
			return super.getDifferenceEstimate(estimate2);
		}
	}

	@Override
	public Estimate<?> getSumEstimate(Estimate<?> estimate2) {
		if (this.isMergeableEstimate(estimate2)) {
			return add((BootstrapHybridPointEstimate) estimate2);
		} else {
			return super.getSumEstimate(estimate2);
		}
	}

	@Override
	public Estimate<?> getProductEstimate(double scalar) {
		return multiply(scalar);
	}

//	@Override
//	public Estimate<?> collapseEstimate(LinkedHashMap<String, List<String>> desiredIndicesForCollapsing) {
//		VarianceEstimate vpe = getVarianceEstimate();
//		Matrix collapsedPointEstimate = collapseRowVector(vpe.getMean(), desiredIndicesForCollapsing);
//		Matrix collapsedVarMean = collapseSquareMatrix(vpe.varMean, desiredIndicesForCollapsing);
//		Matrix collapsedMeanVar = collapseSquareMatrix(vpe.meanVar, desiredIndicesForCollapsing);
//		Matrix collapsedDesignVarianceOfMeanRealizedY = collapseSquareMatrix(vpe.designVarianceOfMeanRealizedY, desiredIndicesForCollapsing);
//		List<String> newIndexRow = new ArrayList<String>(desiredIndicesForCollapsing.keySet());
//		Collections.sort(newIndexRow);
//		VarianceEstimate outputEstimate = new VarianceEstimate(vpe.numberOfRealizations,
//				collapsedPointEstimate,
//				collapsedVarMean,
//				collapsedMeanVar,
//				collapsedDesignVarianceOfMeanRealizedY,
//				newIndexRow);
//		return outputEstimate;
//	}

	
}
