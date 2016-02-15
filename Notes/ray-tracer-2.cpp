#define MAX_RAY_DEPTH 3

color Trace(const Ray &ray, int depth) {

	Object *object = nullptr;
	float minDist = INFINITY;
	Point pHit;
	Normal nHit;

	for (int k = 0; k < objects.size(); k++) {
		if (intersect(objects[k], ray, &pHit, &nHit)) {
			float distance = distance(ray.origin, pHit);
			if (distance < minDist) {
				object = objects[i];
				minDist = distance;
			}
		}
	}

	if (object == nullptr)
		return 0;

	if (object->isGlass && depth < MAX_RAY_DEPTH) {
		Ray reflectionRay = computeReflectionRay(ray.direction, nHit);
		color reflectionColor = Trace(reflectionRay, depth+1);
		Ray refractionRay = computeRefractionRay(object->indexOfRefraction, ray.direction, nHit);
		color refractionColor = Trace(refractionColor, depth+1);
		float Kr, Kt;
		fresnel(object->indexOfRefraction, nHit, ray.direction, &Kr, &Kt);
		return reflectionColor * Kr + refractionColor * (1-Kr);
	}

	Ray shadowRay;
	shadowRay.direction = lightPosition - pHit;
	bool inShadow = false;
	for (int k = 0; k < objects.size(); k++) {
		if (intersect(objects[k], shadowRay))
			return 0;
	}
	
	return object->color * light.brightness;
}

for (int j = 0; j < imageHeight; j++) {
	for (int i = 0; i < imageWidth; i++) {
		Ray primRay;
		computePrimRay(i, j, &primRay);
		pixels[i,j] = Trace(primRay, 0);
	}
}