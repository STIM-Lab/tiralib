#pragma once

#include <Windows.h>
#include "tira/image.h"

#include <sstream>

#include "thorcam/tl_camera_sdk_load.h"
#include "thorcam/tl_camera_sdk.h"
#include "thorcam/tl_mono_to_color_processing_load.h"
#include "thorcam/tl_mono_to_color_processing.h"

void frame_available_callback_global(void* sender, unsigned short* image_buffer, int frame_count, unsigned char* metadata, int metadata_size_in_bytes, void* context);



class ThorlabsCamera {

	int _bit_depth;
	long long _exposure;
	int _gain_range[2];
	double _gain_dB = 6.0;
	int _image_size[2];
	void* _camera_handle = 0;
	bool _sdk_open = false;
	bool _mono2color_sdk_open = false;
	bool _dll_init = false;
	enum TL_CAMERA_SENSOR_TYPE _camera_sensor_type;
	enum TL_COLOR_FILTER_ARRAY_PHASE _color_filter_array_phase;
	float _color_correction_matrix[9];
	float _default_white_balance_matrix[9];
	
	unsigned short* _callback_buffer = 0;


	HANDLE _frame_acquired_event = 0;
	volatile int _is_first_frame_finished = 0;
	void* _mono_to_color_processor_handle = 0;

	void allocate_buffers() {
		_callback_buffer = (unsigned short*)malloc(sizeof(unsigned short) * _image_size[0] * _image_size[1]);

		// color image size will be 3x the size of a mono image
		//_buffer = (unsigned short*)malloc(sizeof(unsigned short) * _image_size[0] * _image_size[1] * 3);
		//_image = tira::image<short>(_image_size[0], _image_size[1], 3);
	}



public:

	void setExposure(double milliseconds) {
		_exposure = (long long)(milliseconds * 1000);
		tl_camera_set_exposure_time(_camera_handle, _exposure);
	}

	double getExposure() {
		return (double)_exposure / 1000.0;
	}

	ThorlabsCamera() {
		_image_size[0] = 256;
		_image_size[1] = 256;
		_exposure = 1000000;						// 1000 ms
	}

	void Init() {
		if (!_dll_init) {
			if (tl_camera_sdk_dll_initialize()) {
				std::cout << "ERROR: Failed to initialize ThorCam DLL" << std::endl;
				return;
			}
			_dll_init = true;
		}

		if (!_sdk_open) {
			if (tl_camera_open_sdk()) {
				std::cout << "ERROR: Failed to open ThorCam SDK" << std::endl;
				return;
			}
			_sdk_open = true;
		}

		char camera_ids[1024];
		camera_ids[0] = 0;

		// Discover cameras.
		if (tl_camera_discover_available_cameras(camera_ids, 1024)) {
			std::cout << "ERROR: Failed to discover cameras - " << tl_camera_get_last_error() << std::endl;
			return;
		}
		if (!strlen(camera_ids)) {
			std::cout << "No cameras found." << std::endl;
			return;
		}
		std::cout << "Camera IDs: " << camera_ids << std::endl;

		// Camera IDs are separated by spaces.
		char* p_space = strchr(camera_ids, ' ');
		if (p_space)
			*p_space = '\0'; // isolate the first detected camera
		char first_camera[256];
		strcpy_s(first_camera, 256, camera_ids);
		std::cout << "Opening camera " << first_camera << "...";
		// Connect to the camera (get a handle to it).
		if (tl_camera_open_camera(first_camera, &_camera_handle)) {
			std::cout << std::endl << "ERROR: Failed to discover cameras - " << tl_camera_get_last_error() << std::endl;
			return;
		}

		std::cout << "handle = 0x" << _camera_handle << std::endl;

		if (tl_mono_to_color_processing_initialize()) {
			std::cout << "ERROR: Failed to open mono/color SDK" << std::endl;
			return;
		}
		_mono2color_sdk_open = true;

		tl_camera_get_camera_sensor_type(_camera_handle, &_camera_sensor_type);
		tl_camera_get_color_filter_array_phase(_camera_handle, &_color_filter_array_phase);
		tl_camera_get_color_correction_matrix(_camera_handle, _color_correction_matrix);
		tl_camera_get_default_white_balance_matrix(_camera_handle, _default_white_balance_matrix);
		tl_camera_get_bit_depth(_camera_handle, &_bit_depth);
		tl_mono_to_color_create_mono_to_color_processor(
			_camera_sensor_type,
			_color_filter_array_phase,
			_color_correction_matrix,
			_default_white_balance_matrix,
			_bit_depth,
			&_mono_to_color_processor_handle);

		tl_camera_set_camera_connect_callback(callback_connect, 0);
		tl_camera_set_camera_disconnect_callback(callback_disconnect, 0);


		tl_camera_set_exposure_time(_camera_handle, _exposure);
		tl_camera_get_gain_range(_camera_handle, &_gain_range[0], &_gain_range[1]);
		if (_gain_range[1] > 0) {																// if this camera supports gain, set it to 6.0 decibels

			int gain_index;
			tl_camera_convert_decibels_to_gain(_camera_handle, _gain_dB, &gain_index);		// set the default gain
			tl_camera_set_gain(_camera_handle, gain_index);
		}

		tl_camera_set_frames_per_trigger_zero_for_unlimited(_camera_handle, 1);				// collect one frame per trigger
		tl_camera_set_frame_available_callback(_camera_handle, frame_available_callback_global, this);		// set the camera callback
		//tl_camera_arm(_camera_handle, 1);														// buffer one frame
		tl_camera_get_image_width(_camera_handle, &_image_size[0]);							// get the current camera image width and height
		tl_camera_get_image_height(_camera_handle, &_image_size[1]);

		allocate_buffers();																		// allocate space for the image buffers
	}

	void reInit() {
		//Destroy();
		//Init();
	}

	void Destroy() {
		if (_camera_handle) {
			if (tl_camera_close_camera(_camera_handle)) {
				std::cout << "ERROR: Failed to close camera - " << tl_camera_get_last_error() << std::endl;
			}
			_camera_handle = 0;
		}
		if (_sdk_open) {
			if (tl_camera_close_sdk()) {
				std::cout << "ERROR: Failed to close ThorCam SDK" << std::endl;
			}
			_sdk_open = false;
		}
		if (_dll_init) {
			if (tl_camera_sdk_dll_terminate()) {
				std::cout << "ERROR: Failed to close ThorCam DLL" << std::endl;
			}
			_dll_init = false;
		}
		if (_camera_handle) {
			tl_camera_disarm(_camera_handle);
		}
	}

	void Snap() {
		_is_first_frame_finished = 0;
		tl_camera_arm(_camera_handle, 1);														// buffer one frame
		if (_camera_handle) {
			tl_camera_issue_software_trigger(_camera_handle);

			for (;;) {
				WaitForSingleObject(_frame_acquired_event, INFINITE);
				if (_is_first_frame_finished) break;
			}
		}
		tl_camera_disarm(_camera_handle);
	}

	static void callback_connect(char* camera_serial_number, enum TL_CAMERA_USB_PORT_TYPE usb_bus_speed, void* context) {
		std::cout << "camera " << camera_serial_number << " connected with bus speed = " << usb_bus_speed << std::endl;
	}

	static void callback_disconnect(char* camera_serial_number, void* context) {
		std::cout << "camera " << camera_serial_number << " disconnected!" << std::endl;
	}

	void frame_available_callback(void* sender, unsigned short* image_buffer, int frame_count, unsigned char* metadata, int metadata_size_in_bytes) {
		if (_is_first_frame_finished)
			return;

		std::cout << "image buffer = 0x" << image_buffer << std::endl;
		std::cout << "frame_count = " << frame_count << std::endl;
		std::cout << "meta data buffer = 0x" << metadata << std::endl;
		std::cout << "metadata size in bytes = " << metadata_size_in_bytes << std::endl;

		_is_first_frame_finished = 1;

		// If you need to save the image data for application specific purposes, this would be the place to copy it into separate buffer.
		if (_callback_buffer)
			memcpy(_callback_buffer, image_buffer, (sizeof(unsigned short) * _image_size[0] * _image_size[1]));

		SetEvent(_frame_acquired_event);

	}

	tira::image<unsigned char> getImage8() {
		tira::image<unsigned char> img8(_image_size[0], _image_size[1], 3);
		if (_camera_handle) {														// If there is an active camera
			tl_mono_to_color_transform_to_24(_mono_to_color_processor_handle,		// transform the buffer into a color image
				_callback_buffer,
				_image_size[0],
				_image_size[1],
				img8.data());
		}
		else {																		// Otherwise fill the texture with a simple color gradient (for debugging)
			for (size_t yi = 0; yi < _image_size[1]; yi++) {
				for (size_t xi = 0; xi < _image_size[0]; xi++) {
					img8(xi, yi, 0) = (unsigned char)((double)xi / (double)_image_size[0] * 256);
					img8(xi, yi, 1) = (unsigned char)((double)yi / (double)_image_size[1] * 256);
					img8(xi, yi, 2) = 128;
				}
			}
		}
		return img8;
	}

	tira::image<unsigned short> getImage16() {
		tira::image<unsigned short> img16(_image_size[0], _image_size[1], 3);
		if (_camera_handle) {
			tl_mono_to_color_transform_to_48(_mono_to_color_processor_handle,
				_callback_buffer,
				_image_size[0],
				_image_size[1],
				img16.data());
		}
		return img16;
	}

};

void frame_available_callback_global(void* sender, unsigned short* image_buffer, int frame_count, unsigned char* metadata, int metadata_size_in_bytes, void* context) {
	((ThorlabsCamera*)context)->frame_available_callback(sender, image_buffer, frame_count, metadata, metadata_size_in_bytes);
}