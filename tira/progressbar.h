#pragma once

namespace tira {
	inline void progressbar(float progress) {
		if (progress > 1.0f) progress = 1.0f;
		int barWidth = 70;

		std::cout << "[";
		int pos = barWidth * progress;
		for (int i = 0; i < barWidth; ++i) {
			if (i < pos) std::cout << "=";
			else if (i == pos) std::cout << ">";
			else std::cout << " ";
		}
		std::cout << "] " << int(progress * 100.0) << " %\r";
		std::cout.flush();
	}
}