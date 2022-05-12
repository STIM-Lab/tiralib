#pragma once

#include "volume.h"
#include "glTexture.h"

namespace tira{

	class glVolume : public volume {

	protected:
		glTexture m_texture;

	};
}