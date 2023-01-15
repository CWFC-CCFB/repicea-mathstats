/*
 * This file is part of the repicea-simulation library.
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
package repicea.simulation.allometrycalculator;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * This calculator provides the basic allometric features, such as the basal area, the mean quadratic diameter, etc.
 * @author Mathieu Fortin - October 2011
 */
public class AllometryCalculator {
	
	public final static double NumberOfTreesHaForDominantFeatures = 100d;
	
	
	/**
	 * This method returns the mean quadratic diameter for a collection of trees. 
	 * @param trees a Collection object that should contain instances of LightAllometryCalculableTree 
	 * @return the mean quadratic diameter in cm
	 */
	public double getMeanQuadraticDiameterCm(Collection<? extends LightAllometryCalculableTree> trees) {
		checkCollection(trees);			
		double d2 = 0;
		double sumND2 = 0;
		double Dg = 0;
		double n = 0;
		double totalN = 0;
		for (LightAllometryCalculableTree t : trees) {
			if (t.getDbhCm() > 0) {
				if (t.getNumber() > 0) {
					d2 = t.getSquaredDbhCm();	// square of dbh
					n = t.getNumber();			// number of stems
					sumND2 += n * d2;
					totalN += n;
				}
			}
		}
		if (trees.size() != 0 && totalN > 0) {
			Dg = Math.sqrt (sumND2 / totalN);
		}
		return Dg;
	}
	
	/**	
	 * This method computes the basal area of a collection of trees. 
	 * @param trees a Collection object that should contain instances of LightAllometryCalculableTree 
	 * @return the basal area in m2 (double)
	 */
	public double getBasalAreaM2(Collection<? extends LightAllometryCalculableTree> trees) {

		checkCollection(trees);

		double basalArea = 0;
		for (LightAllometryCalculableTree t : trees) {
			if (t.getDbhCm() > 0) {
				if (t.getNumber() > 0) {
					basalArea += t.getStemBasalAreaM2() * t.getNumber();
				}
			}
		}
		return basalArea;
	}

	/**
	 * This method returns the number of trees for a collection of trees.
	 * @param trees a Collection object that should contain instances of LightAllometryCalculableTree 
	 * @return the number of trees (double)
	 */
	public double getNumberOfTrees(Collection<? extends LightAllometryCalculableTree> trees) {

		checkCollection(trees);

		double numberOfStems = 0;
		for (LightAllometryCalculableTree t : trees) {
			if (t.getDbhCm() > 0) {
				numberOfStems += t.getNumber();
			}
		}
		return numberOfStems;
	}

	/**
	 * This method returns a collection of CalculableTree objects that have a dbh greater than a particular
	 * threshold.
	 * @param trees a Collection object that should contain instances of AllometryCalculableTree 
	 * @param threshold the dbh threshold in cm (double)
	 * @return a Collection object
	 */
	public Collection<LightAllometryCalculableTree> getTreesLarger(Collection<? extends LightAllometryCalculableTree> trees, double threshold) {

		checkCollection(trees);
		
		Collection<LightAllometryCalculableTree> treesLarger = new ArrayList<LightAllometryCalculableTree>();

		for (LightAllometryCalculableTree t : trees) {
			if (t.getDbhCm() > threshold) {
				treesLarger.add(t);
			}
		}
		return treesLarger;
	}
	
	
	/**
	 * This method computes the commercial volume for a collection of commercial trees.
	 * NOTE: The commercial volume is defined as the  volume for trees that are equal to or greater than 
	 * a particular dbh (e.g. 9.1 cm in dbh for Quebec).
	 * @param trees a Collection object that should contains AllometryCalculableTree instances
	 * @return the commercial volume in m3 (double)
	 */
	public double getCommercialVolumeM3(Collection<? extends AllometryCalculableTree> trees) {

		checkCollection(trees);

		double volume = 0;
		for (AllometryCalculableTree t : trees) {
			if (t.getNumber() > 0) {
				volume += t.getCommercialVolumeM3() * t.getNumber();
			}
		}
		return volume;
	}

	
	/**
	 * This method computes the above ground volume for a collection of trees.
	 * @param trees a Collection object that should contains AllometryCalculableTree instances
	 * @return the total volume in m3 (double)
	 */
	public double getTotalVolumeM3(Collection<? extends AllometryCalculableTree> trees) {

		checkCollection(trees);

		double volume = 0;
		for (AllometryCalculableTree t : trees) {
			if (t.getNumber() > 0) {
				volume += t.getTotalVolumeM3() * t.getNumber();		
			}
		}
		return volume;
	}

	
	
	
	/**
	 * This method returns the dominant height for the trees contained in the trees collection. 
	 * @param trees the Collection object that contains instances of AllometryCalculableTree
	 * @param plotAreaHa the area over which the trees were measured in ha (double)
	 * @param weighted true to enable the plot weighting or false otherwise
	 * @return the dominant height in m (double)
	 */
	public double getDominantHeightM(Collection<? extends AllometryCalculableTree> trees,
			double plotAreaHa,
			boolean weighted) {
		
		checkCollection(trees);

		return getDominantFeature(trees, plotAreaHa, true, weighted);
	}
	
	
	/**
	 * This method returns the dominant diameter for the trees contained in the trees collection. 
	 * @param trees the Collection object that contains instances of AllometryCalculableTree
	 * @param plotAreaHa the area over which the trees were measured in ha (double)
	 * @param weighted true to enable the plot weighting or false otherwise
	 * @return the dominant diameter in cm (double)
	 */
	public double getDominantDiameterCM(Collection<? extends AllometryCalculableTree> trees,
			double plotAreaHa,
			boolean weighted) {
		
		checkCollection(trees);

		return getDominantFeature(trees, plotAreaHa, false, weighted);
	}

	
	
	
	/**
	 * This method returns a dominant feature defined by the method parameter.
	 * @param trees collection of AllometryCalculableTree instances
	 * @param plotAreaHa the plot area in ha (double)
	 * @param height true to get the dominant height or false to get the dominant diameter
	 * @param weighted true to weight the dominant feature by the plot weight or false to use a default weight of 1.0
	 * @return the dominant feature (double)
	 */
	private double getDominantFeature(Collection<? extends AllometryCalculableTree> trees,
			double plotAreaHa,
			boolean height,
			boolean weighted) {		
		
		List<AllometryCalculableTree> copyList = new ArrayList<AllometryCalculableTree>();
		copyList.addAll(trees);
		Collections.sort(copyList, new DbhComparator());

		double domFeature = 0;
		double numberRepresentedByThisTree;
		double numberAdd;
		double numberTreesSoFar = 0;
		double areaFactor = 1d / plotAreaHa;
		double weightingFactor;
		double numberOfTreesPerHectareForDominance = NumberOfTreesHaForDominantFeatures;

		if (plotAreaHa < 0.5) {
			numberOfTreesPerHectareForDominance = (NumberOfTreesHaForDominantFeatures * plotAreaHa - 1) / plotAreaHa;
		} 		
		
		while (!copyList.isEmpty() && numberTreesSoFar < numberOfTreesPerHectareForDominance) {
			AllometryCalculableTree tree = copyList.remove(copyList.size() - 1);	// we remove the last element of the list, ie the largest tree
			if (tree.getNumber() > 0) {
				if (weighted) {
					weightingFactor = tree.getPlotWeight();
				} else {
					weightingFactor = 1d;
				}
				
				numberRepresentedByThisTree = tree.getNumber() * areaFactor * weightingFactor;
				
				if (numberTreesSoFar + numberRepresentedByThisTree <= numberOfTreesPerHectareForDominance) {
					numberAdd = numberRepresentedByThisTree;
				} else {
					numberAdd = numberOfTreesPerHectareForDominance - numberTreesSoFar; // add the remaining part
				}
				
				if (height) {
					domFeature += tree.getHeightM() * numberAdd;
				} else {
					domFeature += tree.getSquaredDbhCm() * numberAdd;
				}
				
				numberTreesSoFar += numberAdd;
			}
		}

		if (numberTreesSoFar > 0) {
			domFeature /= numberTreesSoFar;
			if (height) {
				return domFeature;
			} else {
				return Math.sqrt(domFeature);	// dominant mean diameter is actually the quadratic diameter
			}
		} else {
			return 0d;		// there is no tree at all
		}
		
//		double dom = 0;
//		double numberAdd;
//		double numberTreesSoFar = 0;
//		double areaFactor = 10000d / plotArea;
//		double[][] values = new double[2][trees.size()];
//		double weightingFactor;
//		int iter = 0;
//		if (!trees.isEmpty()) {
//			for (AllometryCalculableTree t : trees) {
//				if (t.getNumber() > 0) {
//					if (weighted) {
//						weightingFactor = t.getPlotWeight();
//					} else {
//						weightingFactor = 1d;
//					}
//					if (height) {
//						values[0][iter] = t.getHeightM();
//					} else {
//						values[0][iter] = t.getDbhCm();
//					}
//					values[1][iter] = t.getNumber() * areaFactor * weightingFactor;
//					iter++;
//				}
//			}
//			double denum = (double) 1 / numberOfTreesPerHectareForDominance;
//			int pointer = -1;
//			while (ObjectUtility.isThereAnyElementDifferentFrom(values[0],-1d) && numberTreesSoFar < numberOfTreesPerHectareForDominance) {
//				pointer = ObjectUtility.findMaxInAnArrayOfDouble(values[0]);
//
//				if (numberTreesSoFar + values[1][pointer] <= numberOfTreesPerHectareForDominance) {
//					numberAdd = values[1][pointer];
//				} else {
//					numberAdd = numberOfTreesPerHectareForDominance
//							- numberTreesSoFar; // add the remaining part
//				}
//
//				numberTreesSoFar += numberAdd;
//
//				dom += values[0][pointer] * numberAdd * denum;
//
//				values[0][pointer] = -1d;
//			} 
//		}
//		return dom;
	}

	private void checkCollection(Collection<? extends LightAllometryCalculableTree> trees) {
		if (trees == null) {
			throw new InvalidParameterException("Collection trees is null!");
		}
	}
	
//	/**
//	 * This method sets the number of trees per hectare on which the dominant features are calculated. By default,
//	 * this number is set to 100 trees per hectare.
//	 * @param numberOfTreesHaForDominantFeatures a double 
//	 */
//	public static void setNumberOfTreesHaForDominantFeatures(double numberOfTreesHaForDominantFeatures) {
//		NumberOfTreesHaForDominantFeatures = numberOfTreesHaForDominantFeatures;
//	}
	
}
