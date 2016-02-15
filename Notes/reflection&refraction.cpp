// Reflection and Refraction

color reflectionColor = computeReflectionColor();
color refractionColor = computeRefractionColor();
float reflectionMixValue;
float refractionMixValue;
fresnel(refractiveIndex, normalHit, primaryRayDirection, 
	&reflectionColor, &refractionMixValue);
color glassBallColorAtHit = reflectionMixValue * reflectionColor + (1-reflectionMixValue) * refractionColor;
