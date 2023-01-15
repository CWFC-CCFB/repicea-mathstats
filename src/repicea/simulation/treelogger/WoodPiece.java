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
package repicea.simulation.treelogger;

import java.io.Serializable;
import java.security.InvalidParameterException;
import java.util.HashMap;
import java.util.Map;

/**	
 * A basic wood piece from a logger. All other pieces class should inherit from this class.
 * @author Mathieu Fortin - April 2010
 */
public abstract class WoodPiece implements Serializable {

	public static enum Property {
		@Deprecated
		totalVolume_m3,
		woodVolume_m3,
		barkVolume_m3,
		juvenileWoodVolume_m3,
		totalMass_kg,
		woodMass_kg,
		woodCarbonMass_kg,
		medianDiameter_cm
	}
	
	
	private static final long serialVersionUID = 20100805L;
	
	private int id = -1;
	private final LoggableTree tree;
	private int rank = -1;
	private boolean withBark;
	private boolean withPith;
	private final LogCategory logCategory;
	private final Map<Property, Double> properties;	

	/**
	 * Constructor based on the loggable tree. This constructor checks if the tree is a Numberable instance. If so, the number of
	 * tree is recorded in the expansionFactor variable.
	 * @param logCategory a TreeLogCategory instance
	 * @param tree a LoggableTree instance
	 */
	private WoodPiece(LogCategory logCategory, LoggableTree tree) {
		this.logCategory = logCategory;
		this.tree = tree;
		properties = new HashMap<Property, Double>();
	}

	
	/**
	 * Constructor based on the loggable tree. This constructor checks if the tree is a Numberable instance. If so, the number of
	 * tree is recorded in the expansionFactor variable.
	 * @param logCategory a TreeLogCategory instance
	 * @param tree a LoggableTree instance
	 * @param id an Integer that represents the id of this piece
	 */
	private WoodPiece(LogCategory logCategory, LoggableTree tree, int id) {
		this(logCategory, tree);
		this.id = id;
	}
	
	
	/**
	 * Constructor based on the loggable tree. This constructor checks if the tree is a Numberable instance. If so, the number of
	 * tree is recorded in the expansionFactor variable.
	 * @param logCategory a TreeLogCategory instance
	 * @param tree a LoggableTree instance
	 * @param id an Integer that represents the id of this piece
	 * @param rank an Integer that is the rank of this log in the tree from the stump to the top
	 */
	private WoodPiece(LogCategory logCategory, LoggableTree tree, int id, int rank) {
		this(logCategory, tree, id);
		this.rank = rank;
	}
	
	
	/**
	 * Constructor based on the loggable tree. This constructor checks if the tree is a Numberable instance. If so, the number of
	 * tree is recorded in the expansionFactor variable.
	 * @param logCategory a TreeLogCategory instance
	 * @param id an Integer that represents the id of this piece
	 * @param tree a LoggableTree instance
	 * @param rank an Integer that is the rank of this log in the tree from the stump to the top
	 * @param withBark a boolean that indicates whether or not the piece was calculated with bark
	 * @param withPith a boolean that indicates whether or not the piece was calculated with pith
	 */
	@Deprecated
	protected WoodPiece(LogCategory logCategory, int id, LoggableTree tree, int rank, boolean withBark, boolean withPith) {
		this(logCategory, tree, id, rank);
		this.withBark = withBark;
		this.withPith = withPith;
	}
	
	
	/**
	 * Constructor #2. 
	 * @param logCategory a TreeLogCategory instance
	 * @param tree a LoggableTree instance
	 * @param overbark a boolean
	 * @param volumeOfThisWoodPieceM3 the volume of wood in a single log (NOTE: without any weighting)
	 */
	protected WoodPiece(LogCategory logCategory, LoggableTree tree, boolean overbark, double volumeOfThisWoodPieceM3) {
		this(logCategory, tree);
		this.withPith = true;
		if (overbark) {
			double volumeUnderBark = volumeOfThisWoodPieceM3 / (1d + tree.getBarkProportionOfWoodVolume());
			setProperty(Property.woodVolume_m3, volumeUnderBark);
			setProperty(Property.barkVolume_m3, volumeOfThisWoodPieceM3 - volumeUnderBark);
		} else {
			setProperty(Property.woodVolume_m3, volumeOfThisWoodPieceM3);
			setProperty(Property.barkVolume_m3, volumeOfThisWoodPieceM3 * tree.getBarkProportionOfWoodVolume());
		}
		this.withBark = getBarkVolumeM3() > 0d;
	}
	
	/**
	 * This method returns the expansion factor for the number of stems. That is the number of stems
	 * represented by the LoggableTree instance
	 * @return a double
	 */
	protected double getNumberOfStemsExpansionFactor() {
		return tree.getNumber();
	}

	/**
	 * This method returns the plot weight expansion factor in the case of unequal sampling probability or 
	 * unequal areas represented by each plot. 
	 * @return a double
	 */
	protected double getPlotWeightExpansionFactor() {
		return tree.getPlotWeight();
	}
	
//	/**
//	 * This method sets the overbark volume of this wood piece, without consideration for any expansion factor. For instance, if the
//	 * LoggableTree is a Numberable instance, this volume correspond to the log coming from a single tree.
//	 * @param volumeOfThisWoodPieceM3 the volume (m3) which is a double
//	 */
//	protected void setVolumeM3(double volumeOfThisWoodPieceM3, boolean overbark) {
//		if (overbark) {
//			setProperty(Property.totalVolume_m3, volumeOfThisWoodPieceM3);
//		} else {
//			setProperty(Property.woodVolume_m3, volumeOfThisWoodPieceM3);
//		}
//	}
	
	/**
	 * This method returns the volume (m3) of this wood piece without any expansion factor.
	 * @return a double
	 */
	public double getWoodVolumeM3() {
		Double volumeM3 = getProperty(Property.woodVolume_m3);
		if (volumeM3 == null) {
			throw new InvalidParameterException("The wood volume has not been set!");
		} else {
			return volumeM3;
		}
	}

	/**
	 * This method returns the volume (m3) of this wood piece without any expansion factor.
	 * @return a double
	 */
	public double getBarkVolumeM3() {
		Double volumeM3 = getProperty(Property.barkVolume_m3);
		if (volumeM3 == null) {
			return 0d;
		} else {
			return volumeM3;
		}
	}

	/**
	 * This method returns the volume (m3) of this wood piece including the bark without any expansion factor.
	 * @return a double
	 */
	public double getTotalVolumeM3() {
		return getWoodVolumeM3() + getBarkVolumeM3();
	}


	
	/**
	 * This method returns the value of a given property.
	 * @param property the name of the property (String)
	 * @return the value of the property (Double)
	 */
	public Double getProperty(Property property) {return properties.get(property);}

	
	/**
	 * This method sets a particular property for this wood piece.
	 * @param property a String that is the property name
	 * @param value a Double that is the value of the property
	 */
	public void setProperty(Property property, Double value) {
		properties.put(property, value);
	}
	
	
	/**
	 * This method returns the identification of this piece.
	 * @return an Integer
	 */
	public int getId() {return id;}

	
	/**
	 * This method returns the tree log category associated to this piece.
	 * @return a TreeLogCategory instance
	 */
	public LogCategory getLogCategory() {return logCategory;}
	
	
	
	/**
	 * This method returns the LoggableTree object from which this piece originates.
	 * @return a LoggableTree instance
	 */
	public LoggableTree getTreeFromWhichComesThisPiece() {return this.tree;}
	
	
	/**
	 * This method return the rank of the log in the tree from the stump to the top. If the rank has not been defined,
	 * it returns -1.
	 * @return an Integer
	 */
	public int getRank() {return rank;}

	
	/**
	 * This method returns the weighted volume of this WoodPiece instance in m3. It should include all expansion factors.
	 * @return the volume in m3 (double)
	 */
	public double getWeightedWoodVolumeM3() {
		return getWoodVolumeM3() * getNumberOfStemsExpansionFactor() * getPlotWeightExpansionFactor();
	}

	/**
	 * This method returns the weighted volume of the bark in m3. It should include all expansion factors.
	 * @return the volume in m3 (double)
	 */
	public double getWeightedBarkVolumeM3() {
		return getBarkVolumeM3() * getNumberOfStemsExpansionFactor() * getPlotWeightExpansionFactor();
	}


	/**
	 * This method returns the weighted volume of the wood and the bark in m3. It should include all expansion factors.
	 * @return the volume in m3 (double)
	 */
	public double getWeightedTotalVolumeM3() {
		return getTotalVolumeM3() * getNumberOfStemsExpansionFactor() * getPlotWeightExpansionFactor();
	}
	/**
	 * This method returns true if the piece was calculated with bark or false otherwise.
	 * @return a boolean
	 */
	public boolean isWithBark() {return withBark;}
	
	
	/**
	 * This method returns true if the piece was calculated with pith or false otherwise.
	 * @return a boolean
	 */
	public boolean isWithPith() {return withPith;}
	
}


