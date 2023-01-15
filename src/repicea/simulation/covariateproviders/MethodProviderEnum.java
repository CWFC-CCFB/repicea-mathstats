package repicea.simulation.covariateproviders;

import repicea.util.REpiceaTranslator;
import repicea.util.REpiceaTranslator.TextableEnum;

public class MethodProviderEnum {
	public static enum VariableForEstimation implements TextableEnum{
		N("Stem density", "Densit\u00E9 d'arbres"),
		G("Basal area", "Surface terri\u00E8re"),
		V("Volume", "Volume"),
		B("Above Ground Biomass", "Biomasse a\\u00E9rienne");

		VariableForEstimation(String englishText, String frenchText) {
			setText(englishText, frenchText);
		}
		
		@Override
		public String toString() {
			return REpiceaTranslator.getString(this);
		}

		@Override
		public void setText(String englishText, String frenchText) {
			REpiceaTranslator.setString(this, englishText, frenchText);
		}
	}
}
