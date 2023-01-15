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

import java.util.ArrayList;
import java.util.List;

import repicea.math.Matrix;
import repicea.stats.data.DataSet;
import repicea.stats.data.StatisticalDataException;
import repicea.stats.data.StatisticalDataStructure;
import repicea.stats.model.CompositeLogLikelihoodWithExplanatoryVariables;
import repicea.stats.model.IndividualLogLikelihood;
import repicea.stats.model.glm.Family;
import repicea.stats.model.glm.Family.GLMDistribution;
import repicea.stats.model.glm.GeneralizedLinearModel;
import repicea.stats.model.glm.LinkFunction;
import repicea.stats.model.glm.LinkFunction.Type;

/**
 * A class implementing the generalized linear model with measurement
 * error on one variable.
 * @author Mathieu Fortin - July 2022
 */
public class GLMWithMeasurementError extends GeneralizedLinearModel {

	protected final GLMMeasErrorDefinition measError;
	
	public GLMWithMeasurementError(DataSet dataSet, String modelDefinition, Matrix startingValues, GLMMeasErrorDefinition measError) {
		super(dataSet, GLMDistribution.Bernoulli, Type.CLogLog, modelDefinition, null, startingValues, measError);
		this.measError = measError;
		this.setConvergenceCriterion(1E-6);
	}
	
	public GLMWithMeasurementError(DataSet dataSet, String modelDefinition, GLMMeasErrorDefinition measError) {
		this(dataSet, modelDefinition, null, measError);
	}

	/*
	 * For extended visibility.
	 */
	@Override
	public StatisticalDataStructure getDataStructure() {
		return super.getDataStructure();
	}

	/*
	 * Just for extended visibility
	 */
	protected LinkFunction getLinkFunction() {return super.getLinkFunction();}
	
	IndividualLogLikelihood getIndividualLogLikelihood() {return individualLLK;}
	
	Matrix getMatrixX() {return matrixX;}
	
	Matrix getVectorY() {return y;}

	@Override
	protected void setModelDefinition(String modelDefinition, Object additionalParm) throws StatisticalDataException {
		super.setModelDefinition(modelDefinition, additionalParm);
		((GLMMeasErrorDefinition) additionalParm).validate(this);
	}
	
	@Override
	protected final CompositeLogLikelihoodWithExplanatoryVariables createCompleteLLK(Object addParm) {
		CompositeLogLikelihoodWithExplanatoryVariables cLL = ((GLMMeasErrorDefinition) addParm).createCompositeLikelihoodFromModel(this);
		return cLL != null ? cLL : super.createCompleteLLK(addParm);
	}
	
	@Override
	protected final StatisticalDataStructure createDataStructure(DataSet dataSet, Object addParm) {
		StatisticalDataStructure sds = ((GLMMeasErrorDefinition) addParm).createDataStructureFromDataSet(dataSet);
		return sds != null ? sds : super.createDataStructure(dataSet, addParm);
	}

	@Override
	protected Family createFamily(GLMDistribution d, Type linkFunctionType, Object addParm) {
		LinkFunction lf = ((GLMMeasErrorDefinition) addParm).createLinkFunction(linkFunctionType, this);
		Family f = lf != null ? 
				Family.createFamily(d, lf) :
					super.createFamily(d, linkFunctionType, null);
		return f;
	}

	@Override 
	protected final IndividualLogLikelihood createIndividualLLK(Object addParm) {
		IndividualLogLikelihood llk = ((GLMMeasErrorDefinition) addParm).createIndividualLogLikelihoodFromModel(this);
		return llk != null ? llk : super.createIndividualLLK(addParm);
	}

	@Override
	public List<String> getEffectList() {
		List<String> effects = new ArrayList<String>();
		effects.addAll(super.getEffectList());
		effects.addAll(measError.getAdditionalEffects());
		return effects;
	}
}
