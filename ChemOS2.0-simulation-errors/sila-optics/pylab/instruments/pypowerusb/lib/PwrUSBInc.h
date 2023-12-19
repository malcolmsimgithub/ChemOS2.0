#ifndef _PWRUSBINC_H_
#define _PWRUSBINC_H_

// Example Commandline application code at the end of this file

// Current Version: 1.3 as of date July 22,2011

// Normal
#ifdef PWRUSBDLL_EXPORTS
   #define CLASS_DECLSPEC   extern "C"  __declspec(dllexport) int  
#else
   #define CLASS_DECLSPEC   extern "C"  __declspec(dllimport) int 
#endif

// With __stdcall

//#ifdef PWRUSBDLL_EXPORTS
//   #define CLASS_DECLSPEC   extern "C"  __declspec(dllexport) int __stdcall 
//#else
//   #define CLASS_DECLSPEC   extern "C"  __declspec(dllimport) int __stdcall
//#endif


#define POWER_USB_MAXNUM 4
// Functions for Version 1 and 2 Firmware PowerUSB. 
///////////////////////////////////////////////////////////
// Main Power Strip functions
CLASS_DECLSPEC InitPowerUSB(int *model, char firmware[]);
CLASS_DECLSPEC ClosePowerUSB();
CLASS_DECLSPEC SetCurrentPowerUSB(int count);
CLASS_DECLSPEC CheckStatusPowerUSB();

CLASS_DECLSPEC SetPortPowerUSB(int port1, int port2, int port3);
CLASS_DECLSPEC SetDefaultStatePowerUSB(int state1, int state2, int port3);


// Functions for Version 2 and later Firmware PowerUSB. 3 ports read and write
///////////////////////////////////////////////////////////
// these functions are available starting with Nov 2010 release of PowerUSB
CLASS_DECLSPEC ReadPortStatePowerUSB(int *state1, int *state2, int *state3);	
CLASS_DECLSPEC ReadDefaultPortStatePowerUSB(int *state1, int *state2, int *state3);	
CLASS_DECLSPEC GetFirmwareVersionPowerUSB(char firmware[]);
CLASS_DECLSPEC GetModelPowerUSB();

// Current sensing and Reset Functions
////////////////////////////////////
CLASS_DECLSPEC ReadCurrentPowerUSB(int *current);
CLASS_DECLSPEC ReadCurrentCumPowerUSB(int *currentCum);
CLASS_DECLSPEC ResetCurrentCounterPowerUSB();
CLASS_DECLSPEC SetCurrentSensRatioPowerUSB(int currentRatio);
CLASS_DECLSPEC SetOverloadPowerUSB(int overload);
CLASS_DECLSPEC GetOverloadPowerUSB();
CLASS_DECLSPEC ResetBoard();
CLASS_DECLSPEC SetCurrentOffset();

// Digital IO functions. Available only when Digital IO model of PowerUSB is connected
////////////////////////////////////
CLASS_DECLSPEC SetIODirectionPowerUSB(int direction[]);
CLASS_DECLSPEC SetOutputStatePowerUSB(int output[]);
CLASS_DECLSPEC GetInputStatePowerUSB(int input[]);
CLASS_DECLSPEC GetInputStateBytePowerUSB(unsigned char *input);
CLASS_DECLSPEC GenerateClockPowerUSB(int port, int onTime, int offTime);
CLASS_DECLSPEC GetOutputStatePowerUSB(int output[]);
CLASS_DECLSPEC SetInputTriggerPowerUSB(int input, int fallingSignal, int outlet, int output, int outTime, char cond, char mask);
CLASS_DECLSPEC SetPLCPowerUSB(int state);
CLASS_DECLSPEC GetPLCPowerUSB(int *state);
CLASS_DECLSPEC ClearPLCPowerUSB();

// Watchdog related functions. Available only when Watchdog model of PowerUSB is connected
//////////////////////////////////////////////////////////////////////////////////////////////
CLASS_DECLSPEC StartWatchdogTimerPowerUSB(int HbTimeSec, int numHbMisses, int resetTimeSec);
CLASS_DECLSPEC StopWatchdogTimerPowerUSB();
CLASS_DECLSPEC GetWatchdogStatusPowerUSB();		// returns 0=IDLE(off), 1=Running, 2=Resetting
CLASS_DECLSPEC SendHeartBeatPowerUSB();
CLASS_DECLSPEC PowerCyclePowerUSB(int resetTimeSec);
CLASS_DECLSPEC ShutdownOffOnPowerUSB(int offDelay, int onDelay);

// SMART Model related Functions
////////////////////////////////////////////////////////////////////////////////////////////////
CLASS_DECLSPEC SetOnOffTimesPowerUSB(int port, int wkDayChunks[], int wkEndChunks[]);
CLASS_DECLSPEC GetOnOffTimesPowerUSB(int port, int wkDayChunks[], int wkEndChunks[]);
CLASS_DECLSPEC SetOnOffFrequencyPowerUSB(int port, int offset, int onTime, int offTime);
CLASS_DECLSPEC SetOnOffModePowerUSB(int port, int mode);
CLASS_DECLSPEC SetModePowerUSB(int mode);
CLASS_DECLSPEC GetModePowerUSB();
CLASS_DECLSPEC SetTVLimitPowerUSB(int wkDayLimit, int wkEndLimit);
CLASS_DECLSPEC SetDateAndTimePowerUSB(unsigned char timeStr[]);
CLASS_DECLSPEC DisplayTextPowerUSB(char text[]);
CLASS_DECLSPEC SetPasswordPowerUSB(char pass[]);
CLASS_DECLSPEC SetSmartParamPowerUSB(int volt, int price, int trigCurrent, int offDelay);
CLASS_DECLSPEC GetMasterOnTimesPowerUSB(int times[]);


#endif




/*
// Command Line PowerUSB Application Code Example
/////////////////////////////////////////////////

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include "..\include\PwrUSBInc.h"

int main(int argc, char* argv[])
{
	int p1=0, p2=0, p3=0, pause=0, model, connected=0;
	int in1, in2, in3;
	char firmware[32];

	// Get command line parameter
	if (argc >= 4)
	{
		p1 = atoi(argv[1]);
		p2 = atoi(argv[2]);
		p3 = atoi(argv[3]);
		if (argc > 4 && argv[4][0] == 'p') 
			pause = 1;
	}
	if (argc < 4 || p1 < 0 || p1 > 1 || p2 < 0 || p2 > 1 || p3 < 0 || p3 >1)
	{
		printf("Usage: pwrusbcmd [0,1] [0,1] [0,1] p (0 for off and 1 for on three outlets) (p for pause at the end of command)\r\n");
		printf("Eg. pwrusbcmd 1 1 1(will switch on all 3 outlets)\r\n");
		printf("Eg. pwrusbcmd 1 0 0(will switch on first outlet and switch off second and third outlet)\r\n");
		printf("Enter to continue");
		getchar();
		return 1;
	}

	// Initialize the PowerUSB or else quit
	if (InitPowerUSB(&model, firmware) > 0)
		connected = CheckStatusPowerUSB();
	if (!connected)
	{
		printf("PowerUSB is not connected");
		printf("Enter to continue");
		getchar();
		return -1;
	}

	// Read and Set the outlets states
	printf("PowerUSB Connected. Model: %d. Version:%s\n", model, firmware);

	ReadPortStatePowerUSB(&in1, &in2, &in3);
	printf("Before Setting Port Status: %d %d %d\n", in1, in2, in3);

	SetPortPowerUSB(p1, p2, p3);

	ReadPortStatePowerUSB(&in1, &in2, &in3);
	printf("After Setting Port Status: %d %d %d\n", in1, in2, in3);

	if (pause)
	{
		printf("Enter to continue");
		getchar();
	}
	return 0;
}

*/