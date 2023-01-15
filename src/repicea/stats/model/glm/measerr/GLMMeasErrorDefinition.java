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

import java.util.List;

import repicea.stats.data.DataSet;
import repicea.stats.data.StatisticalDataStructure;
import repicea.stats.model.CompositeLogLikelihoodWithExplanatoryVariables;
import repicea.stats.model.IndividualLogLikelihood;
import repicea.stats.model.glm.LinkFunction;
import repicea.stats.model.glm.LinkFunction.Type;

/**
 * An interface to ensure an error definition in the GLMWithMeasurementError class.
 * @author Mathieu Fortin - July 2022
 */
public interface GLMMeasErrorDefinition {
	
	public static enum MeasurementErrorModel {Classical, Berkson}
	
	
	/**
	 * Return the measurement error model, that is either classical or Berkson.
	 * @return a MeasurementErrorModel enum
	 */
	public MeasurementErrorModel getMeasurementErrorModel();

	/**
	 * Validate the measurement error model with the generalized linear model. <br>
	 * <br>
	 * The implementation of this method should throw InvalidParameterException if
	 * the measurement error model is invalid. 
	 * 
	 * @param glm a GLMWithMeasurementError instance
	 */
	public void validate(GLMWithMeasurementError glm);
	
	/**
	 * Produce the appropriate CompositeLogLikelihoodWithExplanatoryVariables instance given
	 * the measurement error model.<br>
	 * <br>
	 * If this method returns null, then the method in the super class will be used instead.
	 *
	 * @param glm a GLMWithMeasurementError instance
	 * @return a CompositeLogLikelihoodWithExplanatoryVariables instance
	 */
	public CompositeLogLikelihoodWithExplanatoryVariables createCompositeLikelihoodFromModel(GLMWithMeasurementError glm);
	
	
	/**
	 * Produce the appropriate StatisticalDataStructure instance given
	 * the measurement error model.<br>
	 * <br>
	 * If this method returns null, then the method in the super class will be used instead.
	 *
	 * @param dataSet a DataSet instance
	 * @return a StatisticalDataStructure instance
	 */
	public StatisticalDataStructure createDataStructureFromDataSet(DataSet dataSet);
	
	/**
	 * Produce the appropriate IndividualLogLikelihood instance given
	 * the measurement error model.<br>
	 * <br>
	 * If this method returns null, then the method in the super class will be used instead.
	 *
	 * @param glm a GLMWithMeasurementError instance
	 * @return a IndividualLogLikelihood instance
	 */
	public IndividualLogLikelihood createIndividualLogLikelihoodFromModel(GLMWithMeasurementError glm);

	/**
	 * Produce the appropriate LinkFunction instance given
	 * the measurement error model.<br>
	 * <br>
	 * If this method returns null, then the method in the super class will be used instead.
	 *
	 * @param linkFunctionType a Type enum that stands for the link function type
	 * @param glm a GLMWithMeasurementError instance
	 * @return a LinkFunction instance
	 */
	public LinkFunction createLinkFunction(Type linkFunctionType, GLMWithMeasurementError glm);

	/**
	 * Provide the additional labels for the parameters other than those of the model likelihood.
	 * @return a List of strings
	 */
	public List<String> getAdditionalEffects();

	
}
