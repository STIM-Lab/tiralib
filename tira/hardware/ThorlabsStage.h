#pragma once

#include <string>
#include <queue>
#include <fstream>

#include "thorlabs/Thorlabs.MotionControl.Benchtop.StepperMotor.h"

struct MoveCommand {
	int axis;
	double position;
	double velocity;
	double acceleration;
};

class ThorlabsStage {

	std::string m_serial_num;
	bool m_moving[3] = { false, false, false };
	const int m_polling_rate = 200;
	const double dev2real = 0.000002441406247;
	const double real2dev = 1.0 / dev2real;
	std::queue<MoveCommand> m_command_queue;
	std::ofstream m_logstream;
	bool m_logging = false;
	bool m_active[3] = { false, false, false };	// flag whether each stage is active

	std::string timestamp() {
		std::stringstream stream;
		std::time_t t = std::time(nullptr);
		std::tm tm = *std::localtime(&t);
		stream << std::put_time(&tm, "%F %T") << "  ";
		return stream.str();
	}
	std::string number(double d) {
		std::stringstream stream;
		stream << std::fixed << std::setw(8) << std::setprecision(4) << d;
		return stream.str();
	}

	void ProcessCommandQueue() {
		if (m_command_queue.empty()) return;
		if (isMoving()) return;

		int dev_position = -1;
		int dev_velocity = -1;
		int dev_acceleration = -1;
		int axis = m_command_queue.front().axis;                            // get all of the properties of the move command
		double position = m_command_queue.front().position;
		double velocity = m_command_queue.front().velocity;
		double acceleration = m_command_queue.front().acceleration;

		if (m_active[axis - 1]) {                                                                           // if the specified axis is active
			SBC_GetDeviceUnitFromRealValue(m_serial_num.c_str(), axis, position, &dev_position, 0);         // convert each property to device units
			SBC_GetDeviceUnitFromRealValue(m_serial_num.c_str(), axis, velocity, &dev_velocity, 1);
			SBC_GetDeviceUnitFromRealValue(m_serial_num.c_str(), axis, acceleration, &dev_acceleration, 2);

			SBC_SetVelParams(m_serial_num.c_str(), axis, dev_acceleration, dev_velocity);                   // set the velocity and acceleration
			SBC_MoveToPosition(m_serial_num.c_str(), axis, dev_position);                                   // set the move position

			m_moving[axis - 1] = true;
		}
		if (m_logging)
			m_logstream << timestamp() << "COMMAND: Move Axis " << axis << " to " << number(position) << " mm (" << dev_position << ")" << std::endl;
		m_command_queue.pop();                                      //pop the command off of the queue
	}


public:

	void Init() {

		TLI_InitializeSimulations();

		// Build list of connected device
		short build_list_error = TLI_BuildDeviceList();
		if (build_list_error) {
			std::cout << "ERROR (" << build_list_error << ") building device list: TLI_BuildDeviceList()" << std::endl;
			return;
		}

		// get device list size 
		short n = TLI_GetDeviceListSize();
		if (n != 1) {
			std::cout << "ERROR: " << n << " stage devices found (expecting 1)" << std::endl;
			std::cout << "Entering simulation mode (no stages active)..." << std::endl;
		}
		else {
			// get BBD serial numbers
			char serialNos[100];
			TLI_GetDeviceListByTypeExt(serialNos, 100, 70);

			// output list of matching devices
			char* searchContext = nullptr;
			char* p = strtok_s(serialNos, ",", &searchContext);
			m_serial_num = std::string(p);

			//get device details
			TLI_DeviceInfo deviceInfo;
			// get device info from device
			TLI_GetDeviceInfo(m_serial_num.c_str(), &deviceInfo);
			// get strings from device info structure
			char desc[65];
			strncpy_s(desc, deviceInfo.description, 64);
			desc[64] = '\0';
			char serialNo[9];
			strncpy_s(serialNo, deviceInfo.serialNo, 8);
			serialNo[8] = '\0';

			// open device
			if (SBC_Open(m_serial_num.c_str()) == 0)
			{
				if (!SBC_LoadSettings(m_serial_num.c_str(), 1))
					std::cout << "ERROR: Failed to load settings for X - SBC_LoadSettings()" << std::endl;
				if (!SBC_LoadSettings(m_serial_num.c_str(), 2))
					std::cout << "ERROR: Failed to load settings for Y - SBC_LoadSettings()" << std::endl;
				if (!SBC_LoadSettings(m_serial_num.c_str(), 3))
					std::cout << "ERROR: Failed to load settings for Z - SBC_LoadSettings()" << std::endl;


				// start the device polling at 200ms intervals
				SBC_StartPolling(m_serial_num.c_str(), 1, m_polling_rate);
				SBC_StartPolling(m_serial_num.c_str(), 2, m_polling_rate);
				SBC_StartPolling(m_serial_num.c_str(), 3, m_polling_rate);

				// enable device so that it can move
				Enable();
			}
			m_active[0] = true;
			m_active[1] = true;            // Y stage isn't installed yet
			m_active[2] = true;
		}       // end if the stages are found
	}
	void Destroy() {
		// Un-initialize stages
		if (m_active[0] || m_active[1] || m_active[2]) {
			if (SBC_Close(m_serial_num.c_str())) {
				std::cout << "ERROR: unable to close device - SBC_Close()" << std::endl;
			}
		}
		TLI_UninitializeSimulations();
	}

	void UpdateStatus() {
		WORD messageType;
		WORD messageId;
		DWORD messageData;
		if (m_active[0]) {
			if (SBC_GetNextMessage(m_serial_num.c_str(), 1, &messageType, &messageId, &messageData)) {
				if (messageType == 2 && messageId == 1)
					m_moving[0] = false;
			}
		}

		if (m_active[1]) {
			if (SBC_GetNextMessage(m_serial_num.c_str(), 2, &messageType, &messageId, &messageData)) {
				if (messageType == 2 && messageId == 1)
					m_moving[1] = false;
			}
		}

		if (m_active[2]) {
			if (SBC_GetNextMessage(m_serial_num.c_str(), 3, &messageType, &messageId, &messageData)) {
				if (messageType == 2 && messageId == 1)
					m_moving[2] = false;
			}
		}

		ProcessCommandQueue();

	}
	void Enable() {
		SBC_EnableChannel(m_serial_num.c_str(), 1);
		SBC_EnableChannel(m_serial_num.c_str(), 2);
		SBC_EnableChannel(m_serial_num.c_str(), 3);
	}

	void DisableDisable() {
		SBC_DisableChannel(m_serial_num.c_str(), 1);
		SBC_DisableChannel(m_serial_num.c_str(), 2);
		SBC_DisableChannel(m_serial_num.c_str(), 3);
		if (m_logging) {
			m_logstream << "COMMAND Disable Stages" << std::endl;
			StopLog();
		}
	}

	void StartLog(std::string logfile) {
		m_logstream.open(logfile);
		m_logging = true;
	}
	void StopLog() {
		m_logstream.close();
		m_logging = false;
	}

	// stage status functions
	bool isMoving() {
		return (m_moving[0] || m_moving[1] || m_moving[2]);
	}

	bool Logging() {
		return m_logging;
	}

	double getPosition(int axis) {
		if (m_active[axis - 1]) {
			int position = SBC_GetPosition(m_serial_num.c_str(), axis);
			//double result = (double)m_position[0] * dev2real;

			double result;
			short conversion_error =
				SBC_GetRealValueFromDeviceUnit(m_serial_num.c_str(), axis, position, &result, 0);
			if (conversion_error) {
				std::cout << "ERROR (" << conversion_error << ") converting from device units to real units (axis = " << axis << ")" << std::endl;
			}
			return result;
		}
		return 0.0;
	}
	double getVelocity(int axis) {
		if (axis <= 0 || axis > 3) {
			std::cout << "ERROR: Axis out of range" << std::endl;
			exit(1);
		}
		if (m_active[axis - 1]) {
			int dev_acceleration;
			int dev_velocity;
			SBC_GetVelParams(m_serial_num.c_str(), axis, &dev_acceleration, &dev_velocity);

			double result;
			short conversion_error =
				SBC_GetRealValueFromDeviceUnit(m_serial_num.c_str(), axis, dev_velocity, &result, 1);
			if (conversion_error) {
				std::cout << "ERROR (" << conversion_error << ") converting from device units to real units (axis = " << axis << ")" << std::endl;
			}
			return result;
		}
		return 0.0;
	}
	void setVelocity(int axis, double velocity) {
		if (axis <= 0 || axis > 3) {
			std::cout << "ERROR: Axis out of range" << std::endl;
			exit(1);
		}
		if (m_active[axis - 1]) {
			int dev_acceleration;
			int dev_velocity;
			SBC_GetVelParams(m_serial_num.c_str(), axis, &dev_acceleration, &dev_velocity);
			short conversion_error =
				SBC_GetDeviceUnitFromRealValue(m_serial_num.c_str(), axis, velocity, &dev_velocity, 1);
			if (conversion_error) {
				std::cout << "ERROR (" << conversion_error << ") converting from device units to real units (axis = " << axis << ")" << std::endl;
			}

			SBC_SetVelParams(m_serial_num.c_str(), axis, dev_acceleration, dev_velocity);
		}
	}
	double getAcceleration(int axis) {
		if (axis <= 0 || axis > 3) {
			std::cout << "ERROR: Axis out of range" << std::endl;
			exit(1);
		}
		if (m_active[axis - 1]) {
			int dev_acceleration;
			int dev_velocity;
			SBC_GetVelParams(m_serial_num.c_str(), axis, &dev_acceleration, &dev_velocity);

			double result;
			short conversion_error =
				SBC_GetRealValueFromDeviceUnit(m_serial_num.c_str(), axis, dev_acceleration, &result, 2);
			if (conversion_error) {
				std::cout << "ERROR (" << conversion_error << ") converting from device units to real units (axis = " << axis << ")" << std::endl;
			}
			return result;
		}
		return 0.0;
	}

	void setAcceleration(int axis, double acceleration) {
		if (axis <= 0 || axis > 3) {
			std::cout << "ERROR: Axis out of range" << std::endl;
			exit(1);
		}
		if (m_active[axis - 1]) {
			int dev_acceleration;
			int dev_velocity;
			SBC_GetVelParams(m_serial_num.c_str(), axis, &dev_acceleration, &dev_velocity);
			short conversion_error =
				SBC_GetDeviceUnitFromRealValue(m_serial_num.c_str(), axis, acceleration, &dev_acceleration, 2);
			if (conversion_error) {
				std::cout << "ERROR (" << conversion_error << ") converting from device units to real units (axis = " << axis << ")" << std::endl;
			}

			SBC_SetVelParams(m_serial_num.c_str(), axis, dev_acceleration, dev_velocity);
		}
	}
	size_t getCommandQueueSize() {
		return m_command_queue.size();
	}

	void MoveToX(double x) {
		int device = -1;
		if (m_active[0]) {
			SBC_GetDeviceUnitFromRealValue(m_serial_num.c_str(), 1, x, &device, 0);
			SBC_ClearMessageQueue(m_serial_num.c_str(), 1);
			SBC_MoveToPosition(m_serial_num.c_str(), 1, device);
			m_moving[0] = true;
		}
		if (m_logging)
			m_logstream << timestamp() << "COMMAND: Move Axis 1 to " << number(x) << " mm (" << device << ")" << std::endl;
	}
	void MoveToY(double y) {
		int device = -1;
		if (m_active[1]) {
			SBC_GetDeviceUnitFromRealValue(m_serial_num.c_str(), 2, y, &device, 0);
			SBC_ClearMessageQueue(m_serial_num.c_str(), 2);
			SBC_MoveToPosition(m_serial_num.c_str(), 2, device);
			m_moving[1] = true;
		}
		if (m_logging)
			m_logstream << timestamp() << "COMMAND: Move Axis 2 to " << number(y) << " mm (" << device << ")" << std::endl;

	}
	void MoveToZ(double z) {
		int device = -1;
		if (m_active[2]) {
			SBC_GetDeviceUnitFromRealValue(m_serial_num.c_str(), 3, z, &device, 0);
			SBC_ClearMessageQueue(m_serial_num.c_str(), 3);
			SBC_MoveToPosition(m_serial_num.c_str(), 3, device);
			m_moving[2] = true;
		}
		if (m_logging)
			m_logstream << timestamp() << "COMMAND: Move Axis 3 to " << number(z) << " mm (" << device << ")" << std::endl;
	}

	void MoveTo(double x, double y, double z) {
		MoveToX(x);
		MoveToY(y);
		MoveToZ(z);
	}
	void Home(int axis) {
		if (m_active[axis - 1]) {
			SBC_ClearMessageQueue(m_serial_num.c_str(), axis);
			SBC_Home(m_serial_num.c_str(), axis);
		}
		if (m_logging) {
			m_logstream << timestamp() << "COMMAND Home " << axis << std::endl;
		}
	}

	void QueueMove(int axis, double position, double velocity, double acceleration) {
		MoveCommand new_move;
		new_move.axis = axis;
		new_move.position = position;
		new_move.velocity = velocity;
		new_move.acceleration = acceleration;
		m_command_queue.push(new_move);
	}
};