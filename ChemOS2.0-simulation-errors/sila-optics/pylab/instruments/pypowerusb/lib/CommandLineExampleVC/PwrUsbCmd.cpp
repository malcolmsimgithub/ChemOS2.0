// PwrUsbCmd.cpp : Defines the entry point for the console application.
#include "stdafx.h"							// Comment this for Linux
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <memory.h>
#include <conio.h>
#include "..\include\PwrUSBInc.h"

#define MODEL_BASIC 1
#define MODEL_DIGIO 2
#define MODEL_WATCH 3
#define MODEL_SMART 4

int SendHeartbeat(int interval);
int ConvertInputsToMask(char allInputs[], unsigned char *cond, unsigned char *mask);

int main(int argc, char* argv[])
{
	int p1=0, p2=0, p3=0, pause=0, current=0, model, connected=0;
	int in1, in2, in3, outlet[20], numOutlets=0, i, MaxUnits, r;
	int watchdog=0, digIn=0, digOut=0,  watchOpt[4], ioStates[7], watchOptions=0, digOptions=0;
	int singlePort=0, singleOpt[3], singleOptions=0; 
	int trigAction=0, trigParam[8], numTrigParam=0;
	int input, loToHi, pwrOutlet, output;
	char trigCond[6][6];
	unsigned char cond, mask;
	char firmware[32];
	double totPwr;

	memset(ioStates, 0, sizeof(int)*7);					// triggger state of 3 outputs and 4 inputs
	// Parse the command line parameters
	for (i = 1; i < argc && i < 17; i++)
	{
		if (singlePort)
		{
			if (singleOptions < 3)
				singleOpt[i-2] = atoi(argv[i]);
			singleOptions++;
			continue;
		}
		if (watchdog)									// watchdog cannot be used with other parameters. just read parameters and quit loop
		{
			if (watchOptions < 4)
				watchOpt[i-2] = atoi(argv[i]);
			watchOptions++;
			continue;
		}
		if (trigAction)									// trigger action cannot be used with other param. just read param and quit loop
		{
			if (i < 10)
				trigParam[i-2] = atoi(argv[i]);
			else
			if (i < 16)
				strcpy(trigCond[i-10], argv[i]);
			numTrigParam++;
			continue;
		}
		if (digOut && i <= digOut+3)					// digital output: read next 3 parameters
		{
			r = i-digOut;
			switch (r)
			{
			case 1: ioStates[2] = atoi(argv[i]); break;			// parm1 = out1 = array 2;
			case 2: ioStates[1] = atoi(argv[i]); break;			// parm2 = out2 = array 1;
			case 3: ioStates[0] = atoi(argv[i]); break;			// parm3 = out3 = array 0;
			default: break;
			}
			digOptions++;
			continue;
		}

		// Read param and set variables
		if (argv[i][0] == 'p')					// pause					
			pause = 1;
		else
		if (argv[i][0] == 'c')					// current read
			current = 1;
		else
		if (argv[1][0] == 'w')					// watchdog
			watchdog = 1;
		else
		if (argv[i][0] == 'i')					// input digitalIO
			digIn = 1;
		else
		if (argv[i][0] == 'o')					// output digitalIO
			digOut = i;
		else
		if (argv[1][0] == 't')					// trigger digitalIO
			trigAction = 1;
		else
		if (argv[1][0] == 's')					// trigger digitalIO
			singlePort = 1;
		else
		if (argv[i][0] == '0' || argv[i][0] == '1')		// outlet states
		{
			outlet[i-1] = atoi(argv[i]);
			numOutlets++;
		}
	}

	// Display parameter help message when wrong or no parameters are given 
	if (!current && !watchdog && !digOut && !digIn && !trigAction && !singlePort && (!numOutlets || (numOutlets%3)))
	{
		printf("Usage: pwrusbcmd [0,1] [0,1] [0,1] p c (0 for off and 1 for on) \r\n");
		printf("   Up to 15 outlets. Must be multiple of 3 for up to 5 PowerUSBs\r\n");
		printf("   \'p\' for pause at the end of command\r\n"); 
		printf("   \'c\' for power consumed in primary PowerUSB unit\r\n");
		printf("Eg. pwrusbcmd 1 0 0 (will switch on outlet1 and switch off outlet2 & 3)\r\n");
		printf("Eg. pwrusbcmd 1 1 1 c (will siwtch-on all outlets and report power consumption\r\n");
		printf("pwrusbcmd s [1-5] [1-3] [0,1] (Single port on-off. [Unit] [Port] [On-Off])\r\n");
		printf("Eg. pwrusbcmd s 2 2 1 (will switch on the outlet 2 in PowerUSB unit 2) \r\n");
		printf("WATCHDOG: pwrusbcmd w [x] [y] [z] (Will not work with other parameters)\r\n");
		printf("   x=Time between heartbeats(hb), y:Number of hb misses. z:cpu reset time(sec)\r\n");
		printf("   Starts watchdog and start sending hb. Esc to stop monitor and quit\r\n");
		printf("Eg. pwrusbcmd w 60 3 20. (Send hb every 60sec. Reset for 20 sec if 3 hb fail\r\n");
		printf("DIGITAL IO PLC Controller Options\r\n");
		printf("Output: pwrusbcmd o [0,1] [0,1] [0,1].(Sets the 3 Digial Output States)\r\n");
		printf("Input : pwrusbcmd i. (Returns the status of 4 digital inputs. 0 or 1)\r\n");
		printf("Trig  : pwrusbcmd t [0 or 1-4] [0,1] [-3-32000]..(6 time) [2222]..(6 time)\r\n");
		printf("   Sets trigger actions. Cannot use with other params. See manual for details\r\n");
		printf("   First Param: 0-Clear Table, -1-OffPLC, -2-OnPLC, 1-4: Input Port Action\r\n");
		printf("   Action: -3-noAction -2-toggle -1-on 0-off 1-32K-secs on\r\n");
		printf("   Addional Condition [2222] 2-dont care, 1-has to be on, 0-has to be off\r\n");
		printf("Eg. pwrusbcmd o 1 0 1. (Sets outlet 1 & 3 on and outlet 2 to off\r\n");
		printf("Eg. pwrusbdmd i (Will return 1 1 0 0, if outlet 1 & 2 are on and 3 & 4 are off\r\n");
		printf("Enter to continue");
		getchar();
		return 1;
	}
	if (numOutlets%3)										// make sure outlet parameters are in multiples of 3
		numOutlets = (int)(numOutlets/3)*3;

	if ((MaxUnits=InitPowerUSB(&model, firmware)) > 0)		// Initialize the PowerUSB
		connected = CheckStatusPowerUSB();
	if (!connected)											// quit if no powerusb
	{	
		printf("PowerUSB is not connected\r\n");
		printf("Enter to continue");
		getchar();
		return -1;
	}
	if (MaxUnits > 0)										// set to first pwrusb unit
		SetCurrentPowerUSB(0);
	printf("PowerUSB Connected. Number of Devices:%d Model: %d. Version:%s\n", MaxUnits, model, firmware);		

	if (watchdog)											// Watchdog operation. quit after watchdog function
	{
		if (model != MODEL_WATCH)							// confirm watchdog model is connected
		{
			printf("Watchdog unit not connected\r\n");
			getchar();
			return -1;
		}
		if (watchOpt[0] == 0)								// if parameter is 0 stop the watchdog timer
		{
			StopWatchdogTimerPowerUSB();
			printf("Stopped Watchdog\r\n");
			return 0;
		}
		if (watchOptions != 3)								// if not, make sure 3 parameters are given to start watchdog
		{
			printf("Missing Parameters for Watchdog");
			getchar();
			return -1;
		}
		if (watchOpt[0] < 1 || watchOpt[0] > 1000)			// hearbeat interval can be > 0 and less than 1000 sec
		{
			printf("Heartbeat interval should be higher than 0 and less than 1000 seconds\n");
			getchar();
			return -1;
		}
		printf("Starting watchdog monitor and sending heartbeat\r\n");
		printf("Press Esc to stop watchdog and quit the program\r\n");
		StartWatchdogTimerPowerUSB(watchOpt[0], watchOpt[1], watchOpt[2]);			// start watchdog with given parameters
		r = SendHeartbeat(watchOpt[0]);												// start sending heartbeat in a loop
		return r;
	}

	// DigitalIO Options
	if (digIn || digOut || trigAction)							// Make sure digitalIO model is connected
	{
		if (model != 2)		
		{
			printf("Digital IO PLC unit not connected\r\n");
			getchar();
			return -1;
		}
	}
	if (digIn)													// Get the input states and display 
	{
		r = GetInputStatePowerUSB(ioStates);
		printf("Input States:%d %d %d %d\n", ioStates[6], ioStates[5], ioStates[4], ioStates[3]);
	}
	if (digOut)													// Set digital outputs. Make sure 3 output parameters are given
	{
		if (digOptions != 3)
		{
			printf("Missing Parameters for Digital Output\r\n");
			getchar();
			return -1;
		}
		r = SetOutputStatePowerUSB(ioStates);
		printf("Set output states to: %d %d %d\r\n", ioStates[2], ioStates[1], ioStates[0]);
	}
	if (trigAction)												// Set a trigger action for an input change
	{
		if (trigParam[0] == 0)									// If first parameter is 0 clear full table and return
		{
			ClearPLCPowerUSB();
			printf("Clearted the PLC trigger table\r\n");
			return 0;
		}
		if (trigParam[0] < 0)
		{
			if (trigParam[0] == -1)								// If first parameter is -1 Stop the PLC process inside PowerUSB
			{
				SetPLCPowerUSB(0);
				printf("PLC Stopped in PowerUSB\r\n");
			}
			else 
			if (trigParam[0] == -2)								// If first parameter is -2 Start the PLC process inside PowerUSB
			{
				SetPLCPowerUSB(1);
				printf("PLC Started in PowerUSB\r\n");
			}
			return 0;
		}
		if (numTrigParam != 8 && numTrigParam != 14)			// make sure 2+6 or 2+6+6 parameters are given (input, lo->hi, 6 outlet/ouputs, 6 conditions)
		{
			printf("Please make sure to give 8 parameters (or optional 8+6 additional condition parameters)");
			return -1;
		}
		input = trigParam[0]-1;								// input 1-4 is changed to 0-3
		loToHi = trigParam[1];								// 0: hi->lo 1: lo->hi
		for (i = 2; i < 8; i++)
		{
			pwrOutlet = output = 0;
			cond = mask = 0;
			if (trigParam[i] != -3)							// set actions only for parameters that are not -3(-3=no action)
			{
				if (i <= 4)									// outlet/output param is converted to separate outlet and output
					pwrOutlet = i-1;						// outlet will be 1-3
				else
					output = i-4;							// output will be 1-3
				if (numTrigParam == 14)						// get optoinal condition input strings
					ConvertInputsToMask(trigCond[i-2], &cond, &mask);
				SetInputTriggerPowerUSB(input, loToHi, pwrOutlet, output, trigParam[i], cond, mask);		// convert condition to 2 bytes understood by PLC
				printf("Set Input Trigger for input:%d Oulet:%d Output:%d Action:%d\r\n", input+1, pwrOutlet, output, trigParam[i]);
			}
		}
		return 0;											// return after trigger set operation
	}

	if (singlePort)
	{
		if (singleOptions != 3 || singleOpt[0] < 1 || singleOpt[0] > 5 || singleOpt[1] < 1 || singleOpt[1] > 3 || singleOpt[2] < 0 || singleOpt[2] > 1)
		{
			printf("Wrong parameters for single port on-off.\r\nUsage \"pwrusbcmd [Unit Number(1-5)] [Port Number(1-3)] [On-Off(0,1)]\"\r\n");
			getchar();
			return -1;
		}
		if (singleOpt[0] > MaxUnits)
		{
			printf("Unit number %d is not connected\r\n", singleOpt[0]);
			getchar();
			return -1;
		}
		SetCurrentPowerUSB(singleOpt[0]-1);
		ReadPortStatePowerUSB(&in1, &in2, &in3);
		switch (singleOpt[1])
		{
		case 1: in1 = singleOpt[2]; break;
		case 2: in2 = singleOpt[2]; break;
		case 3: in3 = singleOpt[2]; break;
		default: break;
		}
		SetPortPowerUSB(in1, in2, in3);
		return 0;
	}

	// if outlet 0-1s are set, set the outlets in each (up to 5) connected PowerUSBs
	if (numOutlets > 0)
	{
		ReadPortStatePowerUSB(&in1, &in2, &in3);			// Read and display the state of outlets before change
		printf("Before Setting Port Status: %d %d %d\n", in1, in2, in3);
	}
	if (numOutlets >= 3)									// 1st PowerUSB (default)
	{
		SetPortPowerUSB(outlet[0], outlet[1], outlet[2]);	// set the oulets to on or off as per parameters
		ReadPortStatePowerUSB(&in1, &in2, &in3);
		printf("After Setting Port Status (unit1): %d %d %d\n", in1, in2, in3);
	}
	if (numOutlets >= 6 && MaxUnits >= 2)					// 2nd PowerUSB
	{
		printf("Setting outlets for PowerUSB unit2: %d %d %d\n", outlet[3], outlet[4], outlet[5]);
		SetCurrentPowerUSB(1);
		SetPortPowerUSB(outlet[3], outlet[4], outlet[5]);
	}
	if (numOutlets >= 9 && MaxUnits >= 3)					// 3rd PowerUSB
	{
		printf("Setting outlets for PowerUSB unit3: %d %d %d\n", outlet[6], outlet[7], outlet[8]);
		SetCurrentPowerUSB(2);
		SetPortPowerUSB(outlet[6], outlet[7], outlet[8]);
	}
	if (numOutlets >= 12 && MaxUnits >= 4)					// 4th PowerUSB
	{
		printf("Setting outlets for PowerUSB unit4: %d %d %d\n", outlet[9], outlet[10], outlet[11]);
		SetCurrentPowerUSB(3);
		SetPortPowerUSB(outlet[9], outlet[10], outlet[11]);
	}
	if (numOutlets >= 15 && MaxUnits >= 5)					// 5th PowerUSB
	{
		printf("Setting outlets for PowerUSB unit5: %d %d %d\n", outlet[12], outlet[13], outlet[14]);
		SetCurrentPowerUSB(4);
		SetPortPowerUSB(outlet[12], outlet[13], outlet[14]);
	}

	if (current)												// Read and display power consumption information
	{
		int r1, r2, current, currentCum;

		SetCurrentPowerUSB(0);
		r1 = ReadCurrentPowerUSB(&current);						// Present current consumption reading in milli amps
		r2 = ReadCurrentCumPowerUSB(&currentCum);				// Cumulative curretn consumption in amps/minute
		totPwr = (double)currentCum / 60 * 120 / 1000 / 1000;	// convert to hour(/60), power assuming 120V(V*A), /1000 convert to kw, /1000 convert ma to Amp
		if (r1 >= 0 && r2 >= 0)
			printf("Power Consumed now(ma):%d. Total Power Consumed(kwh)(120VAC):%7.3f\n", current, totPwr);
		else
			printf("Error in Current Reading");
	}

	if (pause)													// pause for user to press enter key
	{
		printf("Enter to continue");
		getchar();
	}
	return 0;
}

// Sends Heartbeat to PowerUSB with an inteterval parameter in seconds. Returns 0 upon comletion of loop
int SendHeartbeat(int interval)
{
	char ch;

	interval = (int)(interval*1000 * 0.8);			// Convert to milli seconds. Convert to 80% of monitor time so that hearbeat is sent quicker than expected
	while (1)
	{
		SendHeartBeatPowerUSB();					// Send heartbeat 
		_sleep(interval);							// Sleep for interval
		if (_kbhit())								// check for key in the input buffer
		{
			ch = _getch();
			if (ch == 27)							// escape key pressed
				break;
		}
	}
	StopWatchdogTimerPowerUSB();					// Stop the watchdog timer since we are stopping the heartbeat
	return 0;
}

// Convert the 2222 characters (4 input condition) into 2 bytes understood by the PLC control
int ConvertInputsToMask(char allInputs[], unsigned char *cond, unsigned char *mask)
{
	char ch[4], ch1[4];								// ch is mask. ch1 is cond
	int i;

	for (i = 0; i < 4; i++)							// create mask. only non 2 char will be enabled
	{
		ch[i] = allInputs[i] - '0';
		ch1[i] = allInputs[i] - '0';
		if (ch[i] == 1 || ch[i] == 0)				// only 0 or 1 will be enabled. rest will be 0 for masking
			ch[i] = 1;
		else
			ch[i] = 0;
		if (ch1[i] == 1)							// 2 will be 0 for condition. rest will be same
			ch1[i] = 1;
		else
			ch1[i] = 0;
	}
	*mask = ((ch[0] << 6) & 0x40) | ((ch[1] << 5) & 0x20) |((ch[2] << 4) & 0x10) | ((ch[3] << 3) & 0x08);
	*cond = ((ch1[0] << 6) & 0x40) | ((ch1[1] << 5) & 0x20) |((ch1[2] << 4) & 0x10) | ((ch1[3] << 3) & 0x08);
	return 0;
}

#ifdef LINUX
#include <termios.h>
#include <unistd.h> 

int _kbhit(void)
{  
	struct termios oldt, newt;  
	int ch;  
	int oldf;   
	tcgetattr(STDIN_FILENO, &oldt);  
	newt = oldt;  
	newt.c_lflag &= ~(ICANON | ECHO);  
	tcsetattr(STDIN_FILENO, TCSANOW, &newt);  
	oldf = fcntl(STDIN_FILENO, F_GETFL, 0);  
	fcntl(STDIN_FILENO, F_SETFL, oldf | O_NONBLOCK);   
	ch = getchar();   
	tcsetattr(STDIN_FILENO, TCSANOW, &oldt);  
	fcntl(STDIN_FILENO, F_SETFL, oldf);   
	if(ch != EOF)  
	{    
		ungetc(ch, stdin);    
		return 1;  
	}   
	return 0;
}

int _getch( ) 
{  
	struct termios oldt,newt;  
	int    ch;  

	tcgetattr( STDIN_FILENO, &oldt );  
	newt = oldt;  
	newt.c_lflag &= ~( ICANON | ECHO );  
	tcsetattr( STDIN_FILENO, TCSANOW, &newt );  
	ch = getchar();  
	tcsetattr( STDIN_FILENO, TCSANOW, &oldt );  
	return ch;
}

#include <sys/select.h>
 
int kbhit2(void)
{
	struct timeval tv;
	fd_set read_fd;
 
	tv.tv_sec=0;
	tv.tv_usec=0;
	FD_ZERO(&read_fd);
	FD_SET(0,&read_fd);
 
	if(select(1, &read_fd, NULL, NULL, &tv) == -1)
		return 0;
 
	if(FD_ISSET(0,&read_fd))
		return 1;
	return 0;
}
 

#endif 