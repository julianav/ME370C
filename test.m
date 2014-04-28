water = importPhase('liquidVapor.xml','water');
water_vap = IdealGasMix('gri30.xml');

set(water,'P',1e5,'T',300);
ug1 = intEnergy_mass(water)
hg1 = enthalpy_mass(water)
rhog1 = density(water)