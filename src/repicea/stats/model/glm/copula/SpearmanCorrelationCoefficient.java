/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2012 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.model.glm.copula;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import repicea.math.Matrix;
import repicea.stats.data.DataBlock;
import repicea.stats.data.HierarchicalSpatialDataStructure;
import repicea.stats.data.HierarchicalStatisticalDataStructure;
import repicea.stats.estimates.CorrelationEstimate;
import repicea.stats.model.glm.copula.CopulaLibrary.DistanceCopula;

/**
 * The SpearmanCorrelationCoefficient class makes it possible to compute Spearman's correlation coefficient. The constructor is empty because the method just offers 
 * methods.
 * @author Mathieu Fortin - September 2011
 */
class SpearmanCorrelationCoefficient {

	
	@SuppressWarnings("rawtypes")
	protected static class SpearmanRankingUnit implements Comparable {
		private double value;
		private double rank;
		private double meanRank;
		
		protected SpearmanRankingUnit(double value) {
			this.value = value;
		}

		protected double getRank() {return rank;}
		protected void setRank(double rank) {this.rank = rank;}
		protected void setMeanRank(double meanRank) {this.meanRank = meanRank;}
		protected double getRankDifference() {return rank - meanRank;}
		
		@Override
		public int compareTo(Object arg0) {
			SpearmanRankingUnit unit = (SpearmanRankingUnit) arg0;
			if (this.value < unit.value) {
				return -1;
			} else if (this.value == unit.value) {
				return 0;
			} else {
				return 1;
			}
		}
	}

	protected static class SpearmanCovariance {
		private SpearmanRankingUnit unitA;
		private SpearmanRankingUnit unitB;
		
		
		protected SpearmanCovariance(SpearmanRankingUnit unitA, SpearmanRankingUnit unitB) {
			this.unitA = unitA;
			this.unitB = unitB;
		}

		protected double getCovarianceValueForThisPair() {return unitA.getRankDifference() * unitB.getRankDifference();}
	}

	
	
	private boolean distanceAvailable;
	
	@SuppressWarnings("unchecked")
	protected CorrelationEstimate[] getSpearmanCorrelationCoefficient(FGMCopulaGLModel model) {
		if (model.getCopula() instanceof DistanceCopula) {
			distanceAvailable = true;
		} else {
			distanceAvailable = false;
		}
	
		Matrix residuals = model.getResiduals();
		List<SpearmanRankingUnit> spearmanRankingUnits = new ArrayList<SpearmanRankingUnit>();
		for (int i = 0; i < residuals.m_iRows; i++) {
			spearmanRankingUnits.add(new SpearmanRankingUnit(residuals.getValueAt(i, 0)));
		}
		
		CorrelationEstimate[] output;
		
		List<SpearmanRankingUnit>[] units;
		List<SpearmanCovariance>[] covariances;
		
		if (distanceAvailable) {
			units = new ArrayList[21];
			covariances = new ArrayList[21];
			output = new CorrelationEstimate[21];
		} else {
			units = new ArrayList[1];
			covariances = new ArrayList[1];
			output = new CorrelationEstimate[1];
		}
		
		for (int i = 0; i < units.length; i++) {
			units[i] = new ArrayList<SpearmanRankingUnit>();
			covariances[i] = new ArrayList<SpearmanCovariance>();
		}
 		
 		SpearmanRankingUnit unitA;
 		SpearmanRankingUnit unitB;
 		HierarchicalStatisticalDataStructure data = model.getDataStructure();
 		Map<String, DataBlock> hierarchicalStructure = data.getHierarchicalStructure();
 		int indexA;
 		int indexB;


 		for (String levelID : hierarchicalStructure.keySet()) {
 			List<Integer> observationIndex = hierarchicalStructure.get(levelID).getIndices();	// index of the observation in this group
 			if (observationIndex != null && !observationIndex.isEmpty()) {
 				int nbObs = observationIndex.size();

 				for (int i = 0; i < nbObs; i++) {
 					indexA = observationIndex.get(i);
 					unitA = spearmanRankingUnits.get(indexA);
 					for (int j = i + 1; j < nbObs; j++) {
 						indexB = observationIndex.get(j);
 						unitB = spearmanRankingUnits.get(indexB);
 						if (distanceAvailable) {
 							double distance = ((HierarchicalSpatialDataStructure) data).getDistancesBetweenObservations().get(0).get(indexA).get(indexB);		// FIXME the 0 must be set to something else here MF2021-12-08
 							int roundedDistance = (int) Math.round(distance);
 							if (roundedDistance < units.length) {
 								if (!units[roundedDistance].contains(unitA)) {
 									units[roundedDistance].add(unitA);
 								}
 								if (!units[roundedDistance].contains(unitB)) {
 									units[roundedDistance].add(unitB);
 								}
 								covariances[roundedDistance].add(new SpearmanCovariance(unitA, unitB));
 							} 
 						} else {
 							if (!units[0].contains(unitA)) {
 								units[0].add(unitA);
 							}
 							if (!units[0].contains(unitB)) {
 								units[0].add(unitB);
 							}
 							covariances[0].add(new SpearmanCovariance(unitA, unitB));
 						}
 					}
 				}
 			}
 		}

		CorrelationEstimate r_hat;
		for (int i = 0; i < output.length; i++) {
			r_hat = new CorrelationEstimate(getCorrelationCoefficientForThisBundle(units[i], covariances[i]), units[i].size());
			output[i] = r_hat;
		}
		
		return output;
	}
	
	
	
	private double getCorrelationCoefficientForThisBundle(List<SpearmanRankingUnit> units, List<SpearmanCovariance> covariances) {
		double variance = getSpearmanVariance(units);
		double covariance = getSpearmanCovariance(covariances);
		return covariance / variance;
	}
	
	
	/**
	 * This method sorts and ranks in ascending order the observations.
	 * @param units a Vector of SpearmanRankingUnit instances
	 */
	@SuppressWarnings("unchecked")
	private void sortAndRank(List<SpearmanRankingUnit> units) {
		Collections.sort(units);
		int i = 0;
		int j = 0;
		double sum;
		while (i < units.size()) {
			sum = 0;
			while (j < units.size() && units.get(i).value == units.get(j).value) {
				sum += ++j;
			}
			double rank = sum / (j - i);
			for (int k = i; k < j; k++) {
				units.get(k).setRank(rank);
			}
			i = j;
		}
	}	
	
	
	
	private double getSpearmanCovariance(List<SpearmanCovariance> covariances) {
		double ssq = 0;
		for (SpearmanCovariance covariance : covariances) {
			ssq += covariance.getCovarianceValueForThisPair();
		}
		return ssq / covariances.size();
	}
	
	
	

	private double getSpearmanVariance(List<SpearmanRankingUnit> units) {
		sortAndRank(units);
		double mean = (double) (units.size() + 1) / 2;		// the mean is actually n * n(+1) / (2 * n)
		double ssq = 0;
		double diff;
		for (SpearmanRankingUnit unit : units) {
			unit.setMeanRank(mean);
			diff = unit.getRankDifference();
			ssq += diff * diff;
		}
		return ssq / units.size();
	}
	
}
