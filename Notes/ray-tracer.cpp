// Ray Tracer

// for each pixel in the image
for (int j = 0; j < imageHeight; j++) {
	for (int i = 0; i < imageWidth; i++) {

		// compute primary ray direction
		Ray primRay;
		computeRay(i, j, &primRay);

		// shot prim ray in the scene and search for the closest intersection
		Point pHit;
		Normal nHit;
		float minDistance = INFINITY;
		Object object = nullptr;

		// for each obejct in the buffer
		for (int k = 0; k < objects.size(); k++) {
			if (intersect(objects[k], primRay, &pHit, &nHit)) {
				float distance = Distance(eyePosition, pHit);
				if (distance < minDistance) {
					object = objects[k];
					minDistance = distance;
				}
			}
		}

		// if the prim ray does hit some object
		if (object != nullptr) {

			// compute illumination
			Ray shadowRay;
			shadowRay.direction = lightPosition - pHit;

			// check if the light shots the pixel
			bool inShadow = false;
			for (int k = 0; k < objects.size(); k++) {
				if (intersect(objects[k], shadowRay)) {
					inShadow = true;
					break;
				}
			}
		}

		// shading with respect to if it's in shadow
		if (inShadow)
			pixels[i][j] = 0;
		else
			pixels[i][j] = object->color * light.brightness;
	}
}