/****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "project.h"
#include "demod_fsk.h"
#include "math.h"
#include "data.h"


int comm_data[WIN_LENGTH];
int data_in_bin[2 * DATA_NUM];
unsigned char data_in_hex[DATA_NUM / 2];
float sigx[CFFT_SIZE*2];
float xcorr[CFFT_SIZE];

int status;     //状态字


void main(void)
{
    int counter=0,flag=0,lfm_sp=0,fsk_sp=0,skip_t=0,findit=0;
    int NOdata=0,start_p=0,rp=0,i=0,readn=0;
    float *sigx_TF,*xcorr_FT;
    float decisionvector[M];
    maxStruct maxval, output;

    f28379d_system_init();

    status=0;   //初始状态
    NOdata=0;   findit=0;   start_p=0;  rp=0;
    lfm_sp=0;   fsk_sp=0;   skip_t=0;   flag=0;
    readn=0;
    for(counter=0;counter<CFFT_SIZE;counter++)
    {
        sigx[2*counter]=0;
        sigx[2*counter+1] = 0;
    }

    while(1){
        //    FILE *fp;
        //    fp = fopen("E:\\test2.bin","rb");
        //    if(fp!=NULL){
        ////        fseek(fp,10,SEEK_SET);
        //        fread(comm_data, sizeof(int),1024,fp);
        ////        fscanf(fp,"%d",counter);
        //        fclose(fp);
        //    }
       readn++;
        if(findit==0){
            for(counter = 0; counter < WIN_LENGTH; counter++)
            {
                int ssp = counter+CFFT_SIZE-WIN_LENGTH;
                sigx[2*ssp] = comm_data[counter];
                sigx[2*ssp+1] = 0;
            }
            sigx_TF=sigx;
            cFFT(sigx_TF);
            xCorr(sigx_TF, lfm_fft);
            xcorr_FT = sigx_TF;
            icFFT(xcorr_FT);
            for(counter = 0; counter < CFFT_SIZE; counter++)
            {
                float real,imag;
                int swp;
                real = xcorr_FT[2*counter];
                imag = xcorr_FT[2*counter+1];
                swp = (counter+CFFT_SIZE/2)%CFFT_SIZE;
                xcorr[swp] =sqrt(real*real+imag*imag);
            }
            maxValue(&maxval, xcorr, CFFT_SIZE);
            flag = isPeak(maxval.Val, xcorr, CFFT_SIZE);
            lfm_sp = (maxval.Loc+CFFT_SIZE/2)%CFFT_SIZE-1;
            if(flag==1 && lfm_sp <WIN_LENGTH)
            {
                status=2;       //当前数据有信号
                findit = 1;
                fsk_sp = lfm_sp + LFM_LENGTH + INTERVAL-10;
                skip_t = (fsk_sp - CFFT_SIZE)/WIN_LENGTH;
                start_p = fsk_sp - CFFT_SIZE;
                start_p = start_p - WIN_LENGTH *skip_t;
            }
            else
            {
                status=1;       //当前数据无信号
                int wp = 2*WIN_LENGTH-CFFT_SIZE;
                for(counter=0;counter<CFFT_SIZE-WIN_LENGTH;counter++)
                {
                    sigx[2*counter]=comm_data[counter+wp];
                    sigx[2*counter+1] = 0;
                }
            }
        }//if(findit==0)
        else if(skip_t !=0){
            status=3;       //跳过状态
            findit=1;
            skip_t--;
        }
        else//当先数据帧包含fsk数据
        {
            for(counter=NOdata;counter<DATA_NUM;counter++)
            {
                if(start_p+NSAMP>=WIN_LENGTH){
                    status=5;       //数据不足，读取下一帧
                    rp=WIN_LENGTH-start_p;
                    for(i=0;i<rp;i++)
                        sigx[i]=comm_data[start_p+i];
                    NOdata=counter;
                    start_p=-rp;
                    break;
                }
                else{//数据足够，开始解调
                    status=4;       //开始解调
                    for(i=rp;i<NSAMP;i++)
                        sigx[i]=comm_data[start_p+i];
                    rp=0;
                    start_p+=NSAMP;
                    decisionvector[0] = SquareLawDetection(sigx, F20Sin, F20Cos, NSAMP);
                    decisionvector[1] = SquareLawDetection(sigx, F22Sin, F22Cos, NSAMP);
                    decisionvector[2] = SquareLawDetection(sigx, F24Sin, F24Cos, NSAMP);
                    decisionvector[3] = SquareLawDetection(sigx, F26Sin, F26Cos, NSAMP);
                    maxValue(&output, decisionvector, M);
                    switch(output.Loc)
                    {
                    case 0:
                        data_in_bin[2 * (DATA_NUM - 1 - counter)]      = 0;
                        data_in_bin[2 * (DATA_NUM - 1 - counter) + 1 ] = 0;
                        break;
                    case 1:
                        data_in_bin[2 * (DATA_NUM - 1 - counter)]      = 0;
                        data_in_bin[2 * (DATA_NUM - 1 - counter) + 1 ] = 1;
                        break;
                    case 2:
                        data_in_bin[2 * (DATA_NUM - 1 - counter)]      = 1;
                        data_in_bin[2 * (DATA_NUM - 1 - counter) + 1 ] = 0;
                        break;
                    case 3:
                        data_in_bin[2 * (DATA_NUM - 1 - counter)]      = 1;
                        data_in_bin[2 * (DATA_NUM - 1 - counter) + 1 ] = 1;
                        break;
                    default:
                        status=8;       //错误
                        break;
                    }//switch(output.Loc)
                }//else{//数据足够，开始解调
            }//for(counter=NOdata;counter<DATA_NUM;counter++)
            if(counter == DATA_NUM){
                bin2hex(data_in_hex, data_in_bin, 2*DATA_NUM);
                status=6;       //解调完成
                NOdata=0;   findit=0;   start_p=0;  rp=0;
                lfm_sp=0;   fsk_sp=0;   skip_t=0;   flag=0;
                break;
            }
        }//else//当先数据帧包含fsk数据
    }// while(1)

    // Step 6. IDLE loop. Just sit and loop forever (optional):
    for(;;)
    {

    }
}

void f28379d_system_init(void)
{
    // Step 1. Initialize System Control:
    // PLL, WatchDog, enable Peripheral Clocks
    // This example function is found in the F2837xD_SysCtrl.c file.
    InitSysCtrl();

#ifdef _STANDALONE
#ifdef _FLASH
    // Send boot command to allow the CPU2 application to begin execution
    IPCBootCPU2(C1C2_BROM_BOOTMODE_BOOT_FROM_FLASH);
#else
    // Send boot command to allow the CPU2 application to begin execution
    IPCBootCPU2(C1C2_BROM_BOOTMODE_BOOT_FROM_RAM);
#endif
#endif

    // Step 2. Initialize GPIO:
    // This example function is found in the F2837xD_Gpio.c file and
    // illustrates how to set the GPIO to it's default state.
    InitGpio(); // Skipped for this example

    // Step 3. Clear all interrupts and initialize PIE vector table:
    // Disable CPU interrupts
    DINT;

    // Initialize the PIE control registers to their default state.
    // The default state is all PIE interrupts disabled and flags
    // are cleared.
    // This function is found in the F2837xD_PieCtrl.c file.
    InitPieCtrl();

    // Disable CPU interrupts and clear all CPU interrupt flags:
    IER = 0x0000;
    IFR = 0x0000;

    // Initialize the PIE vector table with pointers to the shell Interrupt
    // Service Routines (ISR).
    // This will populate the entire table, even if the interrupt
    // is not used in this example.  This is useful for debug purposes.
    // The shell ISR routines are found in F2837xD_DefaultIsr.c.
    // This function is found in F2837xD_PieVect.c.
    InitPieVectTable();

    EALLOW;  // This is needed to write to EALLOW protected registers
    PieVectTable.TIMER0_INT = &cpu_timer0_isr;
    //    PieVectTable.SCIA_RX_INT = &sciaRxFifoIsr;
    //    PieVectTable.SCIB_RX_INT = &scibRxFifoIsr;
    //    PieVectTable.SCIC_RX_INT = &scicRxFifoIsr;
    //    PieVectTable.SCID_RX_INT = &scidRxFifoIsr;
    EDIS;   // This is needed to disable write to EALLOW protected registers
    //    PieCtrlRegs.PIECTRL.bit.ENPIE = 1;   // Enable the PIE block
    //    PieCtrlRegs.PIEIER9.bit.INTx1=1;     // PIE Group 9, INT1
    //    PieCtrlRegs.PIEIER9.bit.INTx2=0;     // PIE Group 9, INT2
    PieCtrlRegs.PIEIER1.bit.INTx7 = 1;         // Enable TINT0 in the PIE: Group 1 __interrupt 7
    IER = 0x100 | M_INT1; // Enable CPU INT
    InitCpuTimers();
    ConfigCpuTimer(&CpuTimer0, 200, 1000000);//0.5s

    // Enable global Interrupts and higher priority real-time debug events:
    EINT;  // Enable Global interrupt INTM
    ERTM;  // Enable Global realtime interrupt DBGM


    EALLOW;
    GpioCtrlRegs.GPADIR.bit.GPIO31 = 1;
    GPIO_SetupPinOptions(34, GPIO_OUTPUT, GPIO_PUSHPULL);
    GPIO_SetupPinMux(34, GPIO_MUX_CPU2, 0);
    //TODO Add code to allow configuration of GPADIR from CPU02 using IPC
    EDIS;
    GpioDataRegs.GPADAT.bit.GPIO31 = 1;// turn off LED

    CpuTimer0Regs.TCR.bit.TSS = 0;//start timer0

}

__interrupt void cpu_timer0_isr(void)
{
    //    CpuTimer0.InterruptCount++;

    CpuTimer0.InterruptCount = ( ++CpuTimer0.InterruptCount ) % 1000000;

    if (CpuTimer0.InterruptCount % 2 == 0) {
        GPIO_WritePin(31, 0);
    } else {
        GPIO_WritePin(31, 1);
    }
    //
    // Acknowledge this __interrupt to receive more __interrupts from group 1
    //
    PieCtrlRegs.PIEACK.all = PIEACK_GROUP1 | PIEACK_GROUP9;
}
