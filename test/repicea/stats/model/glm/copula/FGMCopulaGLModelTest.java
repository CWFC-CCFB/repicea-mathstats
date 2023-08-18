package repicea.stats.model.glm.copula;

import static org.junit.Assert.assertEquals;

import java.util.logging.ConsoleHandler;
import java.util.logging.Level;

import org.junit.BeforeClass;
import org.junit.Test;

import repicea.math.ParameterBound;
import repicea.math.optimizer.LikelihoodOptimizer;
import repicea.stats.data.DataSet;
import repicea.stats.data.StatisticalDataException;
import repicea.stats.estimators.MaximumLikelihoodEstimator;
import repicea.stats.model.glm.Family.GLMDistribution;
import repicea.stats.model.glm.GeneralizedLinearModel;
import repicea.stats.model.glm.LinkFunction.Type;
import repicea.stats.model.glm.copula.CopulaLibrary.DistanceLinkFunctionCopulaExpression;
import repicea.util.ObjectUtility;
import repicea.util.REpiceaLogManager;

public class FGMCopulaGLModelTest {

	@BeforeClass
	public static void doThis() {
		Level l = Level.OFF;
		LikelihoodOptimizer.LOGGER_NAME = MaximumLikelihoodEstimator.LOGGER_NAME;
		ConsoleHandler ch = new ConsoleHandler();
		ch.setLevel(l);
		REpiceaLogManager.getLogger(MaximumLikelihoodEstimator.LOGGER_NAME).setLevel(l);
		REpiceaLogManager.getLogger(MaximumLikelihoodEstimator.LOGGER_NAME).addHandler(ch);		
	}

	/**
	 * This test implements a constant parameter copula model.
	 * @throws Exception
	 */
	@Test
    public void TestWithSimpleCopula() throws Exception {
		double expectedCopulaValue = 0.18072198394164396;
		double expectedLlk = -1077.1623243630413;
		String filename = ObjectUtility.getPackagePath(FGMCopulaGLModelTest.class).concat("donneesR_min.csv");
		DataSet dataSet = new DataSet(filename, true);
		
		GeneralizedLinearModel glm = new GeneralizedLinearModel(dataSet, GLMDistribution.Bernoulli, Type.Logit, "coupe ~ diffdhp + marchand:diffdhp + marchand:diffdhp2 +  essence");
		
		try {
			CopulaExpression copula = new CopulaLibrary.SimpleCopulaExpression(0.0, "IDENT");
			FGMCopulaGLModel copulaModel = new FGMCopulaGLModel(glm, copula);
			copulaModel.doEstimation();
			double actual = copula.getValue();
			assertEquals(expectedCopulaValue, actual, 1E-5);
			double actualLlk = copulaModel.getCompleteLogLikelihood().getValue();
			assertEquals(expectedLlk, actualLlk, 1E-5);
		} catch (StatisticalDataException e) {
			e.printStackTrace();
			throw e;
		}
		
	}

	/**
	 * This test implements a constant parameter copula model.
	 * @throws Exception
	 */
	@Test
    public void TestWithDistanceCopula() throws Exception {
		double expectedCopulaValue = -0.17037528040263983;
		double expectedLlk = -1062.135150078125;
		String filename = ObjectUtility.getPackagePath(FGMCopulaGLModelTest.class).concat("donneesR_min.csv");
		DataSet dataSet = new DataSet(filename, true);
		
		GeneralizedLinearModel glm = new GeneralizedLinearModel(dataSet, GLMDistribution.Bernoulli, Type.Logit, "coupe ~ diffdhp + marchand:diffdhp + marchand:diffdhp2 +  essence");
		
		try {
			CopulaExpression distanceCopula = new CopulaLibrary.DistanceLinkFunctionCopulaExpression(Type.Log, "IDENT", "X + Y", -.15);		// no intercept in the linear expression
			((DistanceLinkFunctionCopulaExpression) distanceCopula).setBounds(0, new ParameterBound(null, 0d));
			FGMCopulaGLModel copulaModel = new FGMCopulaGLModel(glm, distanceCopula);
			copulaModel.setConvergenceCriterion(1E-8);
			((MaximumLikelihoodEstimator) copulaModel.getEstimator()).gridSearch(copulaModel.getParameters().m_iRows - 1, -.25d, -.15d, .01);
			copulaModel.doEstimation();
			double actual = distanceCopula.getParameters().getValueAt(0, 0);
			assertEquals(expectedCopulaValue, actual, 1E-5);
			double actualLlk = copulaModel.getCompleteLogLikelihood().getValue();
			assertEquals(expectedLlk, actualLlk, 1E-5);
		} catch (StatisticalDataException e) {
			e.printStackTrace();
			throw e;
		}
	}
}
