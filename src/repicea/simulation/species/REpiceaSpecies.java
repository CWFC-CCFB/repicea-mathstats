package repicea.simulation.species;

import repicea.simulation.covariateproviders.treelevel.BarkProportionProvider;
import repicea.simulation.covariateproviders.treelevel.BasicWoodDensityProvider;
import repicea.simulation.covariateproviders.treelevel.SpeciesTypeProvider;
import repicea.util.REpiceaTranslator;
import repicea.util.REpiceaTranslator.TextableEnum;

public class REpiceaSpecies {

	public static enum Species implements TextableEnum, 
										SpeciesTypeProvider, 
										BarkProportionProvider,
										BasicWoodDensityProvider {
		Abies_spp(SpeciesType.ConiferousSpecies, 0.40, 0.118, "Fir", "Sapin"),
		Acer_spp(SpeciesType.BroadleavedSpecies, 0.52, 0.109, "Maple", "Erable"),
		Alnus_spp(SpeciesType.BroadleavedSpecies, 0.45, 0.115, "Alder", "Aulne"),
		Betula_spp(SpeciesType.BroadleavedSpecies, 0.51, 0.110, "Birch", "Bouleau"),
		Carpinus_betulus(SpeciesType.BroadleavedSpecies, 0.63, 0.086, "Hornbeam", "Charme"), // this one is from IPCC guidelines 2003
		Castanea_sativa(SpeciesType.BroadleavedSpecies, 0.48, 0.150, "Chestnut", "Chataignier"),	// this one is from IPCC guidelines 2003
		Fagus_sylvatica(SpeciesType.BroadleavedSpecies, 0.58, 0.060, "European beech", "H\u00EAtre europ\u00E9en"),
		Fraxinus_spp(SpeciesType.BroadleavedSpecies, 0.57, 0.160, "Ash", "Fr\u00EAne"),
		Juglans_spp(SpeciesType.BroadleavedSpecies, 0.53, 0.150, "Wallnut", "Noyer"),	// this one is from IPCC guidelines 2003
		Larix_decidua(SpeciesType.ConiferousSpecies, 0.46, 0.140, "Larch", "M\u00E9l\u00E0ze"),
		Picea_abies(SpeciesType.ConiferousSpecies, 0.40, 0.126, "Norway spruce", "Epinette de Norv\u00E0ge"),
		Picea_sitchensis(SpeciesType.ConiferousSpecies, 0.40, 0.126, "Sitka spruce", "Epinette de Sitka"),
		Pinus_pinaster(SpeciesType.ConiferousSpecies, 0.44, 0.160, "Maritime pine", "Pin maritime"),
		Pinus_radiata(SpeciesType.ConiferousSpecies, 0.38, 0.134, "Monterey pine", "Pin de Monterey"),
		Pinus_strobus(SpeciesType.ConiferousSpecies, 0.32, 0.160, "White pine", "Pin blanc"),
		Pinus_sylvestris(SpeciesType.ConiferousSpecies, 0.42, 0.160, "Scots pine", "Pin sylvestre"),
		Populus_spp(SpeciesType.BroadleavedSpecies, 0.35, 0.184, "Poplar", "Peuplier"),
		Prunus_spp(SpeciesType.ConiferousSpecies, 0.49, 0.092, "Cherry", "Cerisier"),
		Pseudotsuga_menziesii(SpeciesType.ConiferousSpecies, 0.45, 0.173, "Douglas fir", "Sapin Douglas"),
		Quercus_spp(SpeciesType.BroadleavedSpecies, 0.58, 0.191, "Oak", "Ch\u00EAne"),
		Salix_spp(SpeciesType.BroadleavedSpecies, 0.45, 0.160, "Willow", "Saule"),
		Thuja_plicata(SpeciesType.ConiferousSpecies, 0.31, 0.106, "Red cedar", "Thuya g\u00E9ant"), // this one is from IPCC guidelines 2003
		Tilia_spp(SpeciesType.BroadleavedSpecies, 0.43, 0.105, "Basswood", "Tilleul"),
		Tsuga_spp(SpeciesType.ConiferousSpecies, 0.42, 0.162, "Hemlock", "Pruche"); // this one is from IPCC guidelines 2003
		
		
		final SpeciesType speciesType;
		final double basicWoodDensity;
		final double barkProportionOfWoodVolume;
		
		Species(SpeciesType speciesType, 
				double basicWoodDensity, 
				double barkProportionOfWoodVolume, 
				String englishName,
				String frenchName) {
			this.speciesType = speciesType;
			this.basicWoodDensity = basicWoodDensity;
			this.barkProportionOfWoodVolume = barkProportionOfWoodVolume;
			setText(englishName, frenchName);
		};

		@Override
		public void setText(String englishText, String frenchText) {
			REpiceaTranslator.setString(this, englishText, frenchText);
		}
		
		public String getLatinName() {return name().replace("_", " ");}
		
		@Override
		public String toString() {return REpiceaTranslator.getString(this);}

		@Override
		public double getBasicWoodDensity() {return basicWoodDensity;}
		
		@Override
		public double getBarkProportionOfWoodVolume() {return barkProportionOfWoodVolume;}		
		
		@Override
		public SpeciesType getSpeciesType() {return speciesType;}
		
		
	}
	
}
