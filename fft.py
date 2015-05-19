import numpy as np
import matplotlib.pyplot as plt
from random import randint
from myhdl import *
class my_fft():
    def ROM(self, i_clk, i_read_adr, o_data, CONTENT):
        DEPTH = len(CONTENT)
        WIDTH = len(o_data) - 1
        r_dout = Signal(intbv(0, -2**WIDTH, 2**WIDTH))
        @always_comb
        def exits():
            o_data.next = r_dout
        @always(i_clk.posedge)
        def seq():
            r_dout.next = CONTENT[i_read_adr]
        return exits, seq

    def DUAL_PORT_RAM(self, i_clk, i_read_adr, i_write_adr, i_we, i_data, o_data):
        DEPTH = max(2**len(i_read_adr), 2**len(i_write_adr))
        WIDTH = len(i_data) - 1
        mem = [Signal(intbv(0, -2**WIDTH, 2**WIDTH)) for i in range(DEPTH)]
        r_dout = Signal(intbv(0, -2**WIDTH, 2**WIDTH))
        
        LEN = len(i_data)
        W = LEN/2 - 1
        
        wdinre = Signal(intbv(0, -2**W, 2**W))
        wdinim =  Signal(intbv(0, -2**W, 2**W))
        wdoutre =  Signal(intbv(0, -2**W, 2**W))
        wdoutim =  Signal(intbv(0, -2**W, 2**W))
        @always_comb
        def f_in_out():
            wdinre.next = i_data[LEN:LEN/2].signed()
            wdinim.next = i_data[LEN/2:0].signed()
            wdoutre.next = o_data[LEN:LEN/2].signed()
            wdoutim.next = o_data[LEN/2:0].signed()
        
        
        @always_comb
        def exits():
            o_data.next = r_dout
            #o_data.next = mem[i_read_adr]
        @always(i_clk.posedge)
        def seq():
            r_dout.next = mem[i_read_adr]
            if i_we == 1:
                mem[i_write_adr].next = i_data
        return exits, seq#, f_in_out

    def RAM(self, i_clk, i_adr, i_we, i_data, o_data):
        DEPTH = 2**len(i_adr)
        WIDTH = len(i_data) - 1
        mem = [Signal(intbv(0, -2**WIDTH, 2**WIDTH)) for i in range(DEPTH)]
        r_dout = Signal(intbv(0, -2**WIDTH, 2**WIDTH))
        @always_comb
        def exits():
            o_data.next = r_dout
        @always(i_clk.posedge)
        def seq():
            r_dout.next = mem[i_adr]
            if i_we == 1:
                mem[i_adr].next = i_data

        return exits, seq
    def Shift_reg(self, i_clk, i_rst, i_en, i_new_data, i_data, o_active, o_data, o_addr, DEPTH):
        WIDTH = len(i_data) -1
        WINDOW_WIDTH = WIDTH - 1
        
        w_ram_din = Signal(intbv(0, -2**WIDTH, 2**WIDTH))
        w_ram_dout = Signal(intbv(0, -2**WIDTH, 2**WIDTH))
        w_rom_dout = Signal(intbv(0, -2**WIDTH, 2**WIDTH))
        
        
        w_din = Signal(intbv(0, -2**WIDTH, 2**WIDTH))
        w_dout = Signal(intbv(0, -2**WIDTH, 2**WIDTH))
        ###
        #FFT_len = DEPTH
        ###
        n = int(np.log2(DEPTH))
        r_addr = Signal(intbv(0,0,2**n))
        r_we = Signal(bool(0))
        
        W_blackman = np.blackman(DEPTH)
        W_blackman = W_blackman / max(W_blackman)
        #CONTENT = tuple([int(np.around(W_blackman[i]*(2**WINDOW_WIDTH))) for i in range(DEPTH)])
        
        CONTENT = tuple([int(1*(2**WINDOW_WIDTH)) for i in range(DEPTH)])
        print "CONTENT",CONTENT
        
        ram = self.RAM(i_clk, r_addr, r_we, w_ram_din, w_ram_dout)        
        rom = self.ROM(i_clk, r_addr, w_rom_dout, CONTENT)
        #ram = self.RAM(i_clk, i_rst, i_en, r_addr, r_addr, r_we, w_din, w_dout)
        bit_reverse_addr = tuple([int(bin(i, len(o_addr))[::-1],2) for i in range(DEPTH)])
        
        r_shift_addr = Signal(intbv(0,0,2**n))
        r_shift_dout = Signal(intbv(0, -2**WIDTH, 2**WIDTH))
        r_shift_act  = Signal(bool(0))
        @always(i_clk.posedge)
        def seq():
            if i_rst == 1:
                r_shift_addr.next = 0
                r_shift_dout.next = 0
                r_shift_act.next = 0
                r_we.next = 0
                r_addr.next = 0
            elif i_en == 1:
                r_shift_addr.next = r_addr
                r_shift_act.next = r_we
                r_shift_dout.next = w_ram_din
#                if r_addr == 0:
#                    r_shift_dout.next = i_data
#                else:
#                    r_shift_dout.next = w_ram_dout
#                
                if i_new_data:
                    r_we.next = 1
                    r_addr.next = 0
                else:
                    if r_addr != r_addr.max - 1:
                        r_addr.next = r_addr + True
                        r_we.next = 1
                    else:
                        r_we.next = 0

        MULT_WIDTH = 2*len(i_data) - 1
        w_mult = Signal(intbv(0,-2**MULT_WIDTH, 2**MULT_WIDTH))
        @always_comb
        def f_mult():
            w_mult.next = (r_shift_dout*w_rom_dout)# >> (WIDTH-1)
        @always_comb
        def comb():
            if r_addr == 0:
                w_ram_din.next = i_data
            else:
                w_ram_din.next = w_ram_dout
        @always_comb
        def f_dout():
            o_data.next = w_mult[len(i_data)+WINDOW_WIDTH:WINDOW_WIDTH].signed()#r_shift_dout
        @always_comb
        def f_active():
            o_active.next = r_shift_act
        @always_comb
        def f_addr():
            o_addr.next = bit_reverse_addr[r_shift_addr]
        return ram, rom, seq, comb, f_dout, f_active, f_addr, f_mult

    def FSM_UNION4(self, i_clk, i_en, o_rd_addr, i_data, o_we, o_data, o_wr_addr,  FFT_SIZE, ORDER, i_log_val):
        ZEROS = Signal(intbv(0,0,2))
        LEN = len(i_data)
        S_LEN = LEN - 1
        MY_LEN = LEN/2
        max_FFT_SIZE = 256
        r_rd_addr = [Signal(intbv(0,0,max_FFT_SIZE)) for i in range(2)]
        r_wr_addr = Signal(intbv(0,0,max_FFT_SIZE))
        r_we = Signal(bool(0))
        L = len(i_data) - 1
        #r_data
        r_i_data = [Signal(intbv(0,-2**L, 2**L)) for i in range(2)]
        r_out_data = [Signal(intbv(0,0, 2**LEN)) for i in range(2)]
        #w_bit_reverse_rd_addr = [Signal(bool(0)) for i in range(int(np.log2(FFT_SIZE)))]
        
        w_bit_reverse_rd_addr = Signal(intbv(0,0, max_FFT_SIZE))
        w_bit_reverse_wr_addr = Signal(intbv(0,0, max_FFT_SIZE))
        #w_bit_reverse_w_addr = [Signal(bool(0)) for i in range(int(np.log2(FFT_SIZE)))]
        w_bit_reverse_w_addr = Signal(intbv(0,0, max_FFT_SIZE/2))
        @always_comb
        def f_addr():
            o_rd_addr.next = w_bit_reverse_rd_addr[len(o_rd_addr):]
        @always_comb
        def f_we():
            o_we.next = r_we
        @always_comb
        def f_wr_addr():
            o_wr_addr.next = w_bit_reverse_wr_addr[len(o_wr_addr):]



        @always_comb
        def f_bitrev_w_addr():
            if   i_log_val == 1:
                #w_bit_reverse_w_addr.next = concat(ZEROS[0], ZEROS[0], ZEROS[0], ZEROS[0], ZEROS[0], ZEROS[0], ZEROS[0])
                w_bit_reverse_w_addr.next = concat(False, False, False, False, False, False, False)
            elif i_log_val == 2:
                w_bit_reverse_w_addr.next = concat(r_rd_addr[1][1], False, False, False, False, False,False)
            elif i_log_val == 3:
                w_bit_reverse_w_addr.next = concat(r_rd_addr[1][1],r_rd_addr[1][2], False,False, False, False, False)
            elif i_log_val == 4:
                w_bit_reverse_w_addr.next = concat(r_rd_addr[1][1],r_rd_addr[1][2], r_rd_addr[1][3], False, False, False, False)
            elif i_log_val == 5:
                w_bit_reverse_w_addr.next = concat(r_rd_addr[1][1],r_rd_addr[1][2], r_rd_addr[1][3], r_rd_addr[1][4], False, False, False)
            elif i_log_val == 6:
                w_bit_reverse_w_addr.next = concat(r_rd_addr[1][1],r_rd_addr[1][2], r_rd_addr[1][3], r_rd_addr[1][4], r_rd_addr[1][5], False, False)
            elif i_log_val == 7:
                w_bit_reverse_w_addr.next = concat(r_rd_addr[1][1],r_rd_addr[1][2], r_rd_addr[1][3], r_rd_addr[1][4], r_rd_addr[1][5], r_rd_addr[1][6], False)
            elif i_log_val == 8:
                w_bit_reverse_w_addr.next = concat(r_rd_addr[1][1],r_rd_addr[1][2], r_rd_addr[1][3], r_rd_addr[1][4], r_rd_addr[1][5], r_rd_addr[1][6], r_rd_addr[1][7])
            else:
                w_bit_reverse_w_addr.next = 0

        @always_comb
        def f_bitrev_data_rd_addr():
            if   i_log_val == 1:
                w_bit_reverse_rd_addr.next = concat(r_rd_addr[0][7], r_rd_addr[0][6], r_rd_addr[0][5], r_rd_addr[0][4], r_rd_addr[0][3],r_rd_addr[0][2], r_rd_addr[0][1],r_rd_addr[0][0])
            elif i_log_val == 2:
                w_bit_reverse_rd_addr.next = concat(r_rd_addr[0][7], r_rd_addr[0][6], r_rd_addr[0][5], r_rd_addr[0][4], r_rd_addr[0][3],r_rd_addr[0][2], r_rd_addr[0][0],r_rd_addr[0][1])
            elif i_log_val == 3:
                w_bit_reverse_rd_addr.next = concat(r_rd_addr[0][7], r_rd_addr[0][6], r_rd_addr[0][5], r_rd_addr[0][4], r_rd_addr[0][3],r_rd_addr[0][0], r_rd_addr[0][1],r_rd_addr[0][2])
            elif i_log_val == 4:
                w_bit_reverse_rd_addr.next = concat(r_rd_addr[0][7], r_rd_addr[0][6], r_rd_addr[0][5], r_rd_addr[0][4], r_rd_addr[0][0],r_rd_addr[0][1], r_rd_addr[0][2],r_rd_addr[0][3])
            elif i_log_val == 5:
                w_bit_reverse_rd_addr.next = concat(r_rd_addr[0][7], r_rd_addr[0][6], r_rd_addr[0][5], r_rd_addr[0][0], r_rd_addr[0][1],r_rd_addr[0][2], r_rd_addr[0][3],r_rd_addr[0][4])
            elif i_log_val == 6:
                w_bit_reverse_rd_addr.next = concat(r_rd_addr[0][7], r_rd_addr[0][6], r_rd_addr[0][0], r_rd_addr[0][1], r_rd_addr[0][2],r_rd_addr[0][3], r_rd_addr[0][4],r_rd_addr[0][5])
            elif i_log_val == 7:
                w_bit_reverse_rd_addr.next = concat(r_rd_addr[0][7], r_rd_addr[0][0], r_rd_addr[0][1], r_rd_addr[0][2], r_rd_addr[0][3],r_rd_addr[0][4], r_rd_addr[0][5],r_rd_addr[0][6])
            elif i_log_val == 8:
                w_bit_reverse_rd_addr.next = concat(r_rd_addr[0][0], r_rd_addr[0][1], r_rd_addr[0][2], r_rd_addr[0][3], r_rd_addr[0][4],r_rd_addr[0][5], r_rd_addr[0][6],r_rd_addr[0][7])
            else:
                w_bit_reverse_rd_addr.next = 0

        @always_comb
        def f_bitrev_data_wr_addr():
            if   i_log_val == 1:
                w_bit_reverse_wr_addr.next = concat(r_wr_addr[7], r_wr_addr[6], r_wr_addr[5], r_wr_addr[4], r_wr_addr[3],r_wr_addr[2], r_wr_addr[1],r_wr_addr[0])
            elif i_log_val == 2:
                w_bit_reverse_wr_addr.next = concat(r_wr_addr[7], r_wr_addr[6], r_wr_addr[5], r_wr_addr[4], r_wr_addr[3],r_wr_addr[2], r_wr_addr[0],r_wr_addr[1])
            elif i_log_val == 3:
                w_bit_reverse_wr_addr.next = concat(r_wr_addr[7], r_wr_addr[6], r_wr_addr[5], r_wr_addr[4], r_wr_addr[3],r_wr_addr[0], r_wr_addr[1],r_wr_addr[2])
            elif i_log_val == 4:
                w_bit_reverse_wr_addr.next = concat(r_wr_addr[7], r_wr_addr[6], r_wr_addr[5], r_wr_addr[4], r_wr_addr[0],r_wr_addr[1], r_wr_addr[2],r_wr_addr[3])
            elif i_log_val == 5:
                w_bit_reverse_wr_addr.next = concat(r_wr_addr[7], r_wr_addr[6], r_wr_addr[5], r_wr_addr[0], r_wr_addr[1],r_wr_addr[2], r_wr_addr[3],r_wr_addr[4])
            elif i_log_val == 6:
                w_bit_reverse_wr_addr.next = concat(r_wr_addr[7], r_wr_addr[6], r_wr_addr[0], r_wr_addr[1], r_wr_addr[2],r_wr_addr[3], r_wr_addr[4],r_wr_addr[5])
            elif i_log_val == 7:
                w_bit_reverse_wr_addr.next = concat(r_wr_addr[7], r_wr_addr[0], r_wr_addr[1], r_wr_addr[2], r_wr_addr[3],r_wr_addr[4], r_wr_addr[5],r_wr_addr[6])
            elif i_log_val == 8:
                w_bit_reverse_wr_addr.next = concat(r_wr_addr[0], r_wr_addr[1], r_wr_addr[2], r_wr_addr[3], r_wr_addr[4],r_wr_addr[5], r_wr_addr[6],r_wr_addr[7])
            else:
                w_bit_reverse_wr_addr.next = 0



        W_mod = 2**(MY_LEN-2)
        cos_table = tuple([int(np.around(W_mod*np.cos(2*np.pi*i/(max_FFT_SIZE)))) for i in range(max_FFT_SIZE)])
#        W_cos = tuple(cos_table[i] for i in range(FFT_SIZE/2))
        sin_table = tuple([int(np.around(-W_mod*np.sin(2*np.pi*i/(max_FFT_SIZE)))) for i in range(max_FFT_SIZE)])
#        W_sin = tuple(sin_table[i] for i in range(FFT_SIZE/2))
        W_cos_sin = tuple([int((65536*cos_table[i]) | (sin_table[i] & 0xFFFF)) for i in range(max_FFT_SIZE/2)])
        print '///'
        #print W_cos
        #print W_sin
        print W_cos_sin
        L_W = MY_LEN - 1
        L_MULT = MY_LEN - 1 + 1
        L_MULT_SUM = L_MULT + 1
        w_mult_re = Signal(intbv(0, -2**L_W, 2**L_W))
        w_mult_im = Signal(intbv(0, -2**L_W, 2**L_W))
#        w_sum_re_plus = Signal(intbv(0, -2**L_MULT, 2**L_MULT))
#        w_sum_im_plus = Signal(intbv(0, -2**L_MULT, 2**L_MULT))
#        w_sum_re_minus = Signal(intbv(0, -2**L_MULT, 2**L_MULT))
#        w_sum_im_minus = Signal(intbv(0, -2**L_MULT, 2**L_MULT))
        
        r_Wre_Wim = Signal(intbv(0,-2**S_LEN, 2**S_LEN))
        r_Wre = Signal(intbv(0,-2**L_W, 2**L_W))
        r_Wim = Signal(intbv(0,-2**L_W, 2**L_W))
        r_sum_re_plus = Signal(intbv(0, -2**L_MULT, 2**L_MULT))
        r_sum_im_plus = Signal(intbv(0, -2**L_MULT, 2**L_MULT))
        r_sum_re_minus = Signal(intbv(0, -2**L_MULT, 2**L_MULT))
        r_sum_im_minus = Signal(intbv(0, -2**L_MULT, 2**L_MULT))
        
#        w_Wre = Signal(intbv(0,-2**L_W, 2**L_W))
#        w_Wim = Signal(intbv(0,-2**L_W, 2**L_W))
#
#        @always_comb
#        def f_WRe():
#            w_Wre.next = W_cos[w_bit_reverse_w_addr]
#        @always_comb
#        def f_WIm():
#            w_Wim.next = W_sin[w_bit_reverse_w_addr]
        
        w_Dre = Signal(intbv(0,-2**L_W, 2**L_W))
        w_Dim = Signal(intbv(0,-2**L_W, 2**L_W))
        @always_comb
        def f_DRe():
            w_Dre.next = i_data[LEN:MY_LEN].signed()
        @always_comb
        def f_DIm():
            w_Dim.next = i_data[MY_LEN:0].signed()


        w_Rre = Signal(intbv(0,-2**L_W, 2**L_W))
        w_Rim = Signal(intbv(0,-2**L_W, 2**L_W))
        @always_comb
        def f_RRe():
            w_Rre.next = r_i_data[1][LEN:MY_LEN].signed()
        @always_comb
        def f_RIm():
            w_Rim.next = r_i_data[1][MY_LEN:0].signed()
        
        MULT_L_W = 2*MY_LEN - 1
        r_Dre_Wre = Signal(intbv(0, -2**MULT_L_W, 2**MULT_L_W))
        r_Dre_Wim = Signal(intbv(0, -2**MULT_L_W, 2**MULT_L_W))
        r_Dim_Wim = Signal(intbv(0, -2**MULT_L_W, 2**MULT_L_W))
        r_Dim_Wre = Signal(intbv(0, -2**MULT_L_W, 2**MULT_L_W))

#
#        @instance
#        def f_monitor():
#            while True:
#                yield i_clk.negedge
#                if i_en == 1:
#                    if ORDER == 4:
#                        print ">>>>>>", r_rd_addr[1], '/', w_Dre, w_Dim
##                    if ORDER == 4 and o_wr_addr == 1:
##                        print '>>>>>>>>', o_wr_addr,'/', o_data[16:].signed(), o_data[:16].signed()                        
##                        print w_Dre, w_Wre, w_Dim, w_Wim
##                        print w_Dre_Wre, w_Dim_Wim, w_Dre_Wim, w_Dim_Wre
##                        print w_Dre_Wre_signed, w_Dim_Wim_signed
##                        print w_Dre_Wre[16+14:14].signed(), w_Dim_Wim[16+14:14].signed()
##                        print w_Dre_Wre_signed[16+14:14].signed(), w_Dim_Wim_signed[16+14:14].signed()
##                        
##                        print w_Dre_Wim[16+14:14].signed(), w_Dim_Wre[16+14:14].signed()
##                        print w_Dre_Wim_signed[16+14:14].signed(), w_Dim_Wre_signed[16+14:14].signed()
                        
#######################################################################
#        w_Dre_Wre_signed = Signal(intbv(0,-2**MULT_L_W, 2**MULT_L_W))
#        w_Dre_Wim_signed = Signal(intbv(0,-2**MULT_L_W, 2**MULT_L_W))
#        w_Dim_Wim_signed = Signal(intbv(0,-2**MULT_L_W, 2**MULT_L_W))
#        w_Dim_Wre_signed = Signal(intbv(0,-2**MULT_L_W, 2**MULT_L_W))
#        @always_comb
#        def f_prepare_mult():
#            if r_Dre_Wre[31] == 1:
#                w_Dre_Wre_signed.next = r_Dre_Wre + 8192
#            else:
#                w_Dre_Wre_signed.next = r_Dre_Wre
#            if r_Dim_Wim[31] == 1:
#                w_Dim_Wim_signed.next = r_Dim_Wim + 8192
#            else:
#                w_Dim_Wim_signed.next = r_Dim_Wim
#            if r_Dre_Wim[31] == 1:
#                w_Dre_Wim_signed.next = r_Dre_Wim + 8192
#            else:
#                w_Dre_Wim_signed.next = r_Dre_Wim
#            if r_Dim_Wre[31] == 1:
#                w_Dim_Wre_signed.next = r_Dim_Wre + 8192
#            else:
#                w_Dim_Wre_signed.next = r_Dim_Wre
#        
        @always_comb
        def f_mult():
            w_mult_re.next = (r_Dre_Wre[16+14:14].signed())-(r_Dim_Wim[16+14:14].signed())
            w_mult_im.next = (r_Dre_Wim[16+14:14].signed())+(r_Dim_Wre[16+14:14].signed())
#            w_mult_re.next = (w_Dre_Wre_signed[16+14:14].signed())-(w_Dim_Wim_signed[16+14:14].signed())
#            w_mult_im.next = (w_Dre_Wim_signed[16+14:14].signed())+(w_Dim_Wre_signed[16+14:14].signed())

        @always(i_clk.posedge)
        def f_sum():
            if i_en == 0:
                r_sum_re_plus.next = 0
                r_sum_im_plus.next = 0
                r_sum_re_minus.next = 0
                r_sum_im_minus.next = 0
            else:
                r_sum_re_plus.next = w_Rre + w_mult_re
                r_sum_im_plus.next = w_Rim + w_mult_im
                r_sum_re_minus.next = w_Rre - w_mult_re
                r_sum_im_minus.next = w_Rim - w_mult_im

        @always_comb
        def f_data():
            o_data.next = r_out_data[r_wr_addr[0]].signed()
        
        w_out_re0 = Signal(intbv(0, -2**L_MULT, 2**L_MULT))
        w_out_re1 = Signal(intbv(0, -2**L_MULT, 2**L_MULT))
        w_out_im0 = Signal(intbv(0, -2**L_MULT, 2**L_MULT))
        w_out_im1 = Signal(intbv(0, -2**L_MULT, 2**L_MULT))
        @always_comb
        def f_out_data():
            if r_sum_re_plus[16] == 1:
                w_out_re0.next = r_sum_re_plus + True
            else:
                w_out_re0.next = r_sum_re_plus
            
            if r_sum_im_plus[16] == 1:
                w_out_im0.next = r_sum_im_plus + True
            else:
                w_out_im0.next = r_sum_im_plus
            
            if r_sum_re_minus[16] == 1:
                w_out_re1.next = r_sum_re_minus + True
            else:
                w_out_re1.next = r_sum_re_minus
            
            if r_sum_im_minus[16] == 1:
                w_out_im1.next = r_sum_im_minus + True
            else:
                w_out_im1.next = r_sum_im_minus
            
        
        @always(i_clk.posedge)
        def seq():
            if i_en == 0:
                r_i_data[0].next = 0
                r_i_data[1].next = 0
                r_wr_addr.next = 0
                r_rd_addr[0].next = 0
                r_rd_addr[1].next = 0
                r_out_data[0].next = 0
                r_out_data[1].next = 0
                r_we.next = 0
                r_Wre.next = 0
                r_Wim.next = 0
                r_Dre_Wre.next = 0
                r_Dim_Wim.next = 0
                r_Dre_Wim.next = 0
                r_Dim_Wre.next = 0
                r_Wre_Wim.next = 0
            else:
                r_rd_addr[1].next = r_rd_addr[0]
                if r_rd_addr[0] != FFT_SIZE - 1:
                    r_rd_addr[0].next = r_rd_addr[0] + True

                if r_rd_addr[1] == 3:
                    r_we.next = 1
                if r_wr_addr == FFT_SIZE - 1:
                    r_we.next = 0

                if r_we == 1:
                    r_wr_addr.next = r_wr_addr + True

                if r_rd_addr[1][0] == 0:
                    r_i_data[0].next = i_data
                    #r_Wre.next = W_cos[w_bit_reverse_w_addr]
                    #r_Wim.next = W_sin[w_bit_reverse_w_addr]
                    r_Wre_Wim.next = W_cos_sin[w_bit_reverse_w_addr]
                else:
                    r_Dre_Wre.next = w_Dre*r_Wre_Wim[LEN:MY_LEN].signed()
                    r_Dim_Wim.next = w_Dim*r_Wre_Wim[MY_LEN:0].signed()
                    r_Dre_Wim.next = w_Dre*r_Wre_Wim[MY_LEN:0].signed()
                    r_Dim_Wre.next = w_Dim*r_Wre_Wim[LEN:MY_LEN].signed()
                    r_i_data[1].next = r_i_data[0]
                    r_out_data[0].next = concat(w_out_re0[17:1], w_out_im0[17:1])
                    r_out_data[1].next = concat(w_out_re1[17:1], w_out_im1[17:1])

        inst = instances()
        return inst#f_bitrev_data_addr, f_bitrev_w_addr,f_bitrev_data_wr_addr, f_out_data, f_addr, seq, f_mult, f_sum, f_data, f_DRe, f_DIm, f_RRe, f_RIm, f_we,f_wr_addr#, f_WRe, f_WIm#, f_prepare_mult
        
        
    def FFT_FSM(self, i_clk, i_rst, i_en, i_new_data, i_data, o_dout, o_active_write, o_addr, FFT_SIZE):
        #bit_reverse_addr = tuple([int(bin(i, len(o_addr))[::-1],2) for i in range(FFT_SIZE)])
        
        State = enum('IDLE','IDLE_STOP', 'WRITE', 'BF2', 'UNION', encoding='one_hot')
        r_state = Signal(State.IDLE)
        r_cnt = Signal(modbv(0,0,FFT_SIZE/4))
        
        WIDTH = len(i_data) - 1
        
        o_w_shift_active = Signal(bool(0))
        o_w_shift_data = Signal(intbv(0,-2**WIDTH, 2**WIDTH))
        o_w_shift_addr = Signal(intbv(0,0,FFT_SIZE))
        fsm_input = self.Shift_reg(i_clk, i_rst, i_en, i_new_data, i_data, o_w_shift_active, o_w_shift_data, o_w_shift_addr, FFT_SIZE)
        
        L = 2*len(i_data) - 1
        i_w_muxed_data = Signal(intbv(0,-2**L, 2**L))
        i_w_muxed_rd_addr = Signal(intbv(0,0,FFT_SIZE))
        i_w_muxed_wr_addr = Signal(intbv(0,0,FFT_SIZE))
        i_w_muxed_we = Signal(bool(0))

        
        o_w_bf2_data = Signal(intbv(0,-2**WIDTH, 2**WIDTH))
        o_w_bf2_rd_addr = Signal(modbv(0,0,FFT_SIZE))
        o_w_bf2_wr_addr = Signal(modbv(0,0,FFT_SIZE))
        o_w_bf2_we = Signal(bool(0))
        i_r_dualram_read_addr = [Signal(modbv(0,0,FFT_SIZE)) for i in range(2)]
        o_w_union4_rd_addr =  Signal(intbv(0,0,FFT_SIZE))
        o_w_union4_wr_addr =  Signal(intbv(0,0,FFT_SIZE))
        o_w_union4_data = Signal(intbv(0,-2**L, 2**L))
        o_w_union4_we = Signal(bool(0))

        @always_comb
        def f_muxed_data():
            if r_state == State.WRITE:
                i_w_muxed_data.next = (o_w_shift_data << 16)
            else:
                i_w_muxed_data.next = o_w_union4_data

        @always_comb
        def f_muxed_rd_addr():
            if r_state == State.WRITE:
                i_w_muxed_rd_addr.next = 0
            else:
                i_w_muxed_rd_addr.next = o_w_union4_rd_addr
        
    
        @always_comb
        def f_muxed_wr_addr():
            if r_state == State.WRITE:
                i_w_muxed_wr_addr.next = o_w_shift_addr
            else:
                i_w_muxed_wr_addr.next = o_w_union4_wr_addr
    
        @always_comb
        def f_muxed_we():
            if r_state == State.WRITE:
                i_w_muxed_we.next = o_w_shift_active
            else:
                i_w_muxed_we.next = o_w_union4_we
        
        
        
        o_w_dualram_data = Signal(intbv(0,-2**L,2**L))
        fsm_dualram = self.DUAL_PORT_RAM(i_clk, i_w_muxed_rd_addr, i_w_muxed_wr_addr, i_w_muxed_we, i_w_muxed_data, o_w_dualram_data)
        r_en_un4 = Signal(bool(0))
        r_log_val = Signal(intbv(0,0,1+int(np.log2(256))))
        
        @always(i_clk.posedge)
        def seq():
            if i_rst == 1:
                r_state.next = State.IDLE
                r_cnt.next = 0
                r_log_val.next = 0
                r_en_un4.next = 0
            else:
                if r_state == State.IDLE:                
                    if i_new_data:
                        r_cnt.next = r_cnt + True
                        if r_cnt == r_cnt.max - 1:
                            r_state.next = State.IDLE_STOP
                elif r_state == State.IDLE_STOP:
                    r_cnt.next = 0
                    r_state.next = State.WRITE
       
                elif r_state == State.WRITE:
                    if i_w_muxed_wr_addr == i_w_muxed_wr_addr.max - 1:
                        r_state.next = State.UNION
                        r_log_val.next = 1
                    #print [[int(mem[i][16:].signed()), int(mem[i][32:16].signed())] for i in range(len(mem))]
                elif r_state == State.UNION:
                    r_en_un4.next = 1
                    if i_w_muxed_wr_addr == i_w_muxed_wr_addr.max - 1:
                        r_en_un4.next = 0
                        if r_log_val == 8:#int(np.log2(FFT_SIZE)):
                            r_state.next = State.IDLE
                            r_log_val.next = 0
                        else:
                            r_state.next = State.UNION
                            r_log_val.next = r_log_val + True
                else:
                    r_state.next = State.IDLE
                    r_cnt.next = 0
                    r_log_val.next = 0
                    r_en_un4.next = 0
                        
#                elif r_state == State.UNION4:
#                    r_log_val.next = 2
#                    r_en_un4.next = 1
#                    if i_w_muxed_wr_addr == i_w_muxed_wr_addr.max - 1:
#                        r_en_un4.next = 0
#                        r_state.next = State.UNION8
#                    #r_state.next = State.IDLE
#                elif r_state == State.UNION8:
#                    r_log_val.next = 3
#                    r_en_un4.next = 1
#                    if i_w_muxed_wr_addr == i_w_muxed_wr_addr.max - 1:
#                        r_en_un4.next = 0
#                        r_state.next = State.UNION16
#                elif r_state == State.UNION16:
#                    r_log_val.next = 4
#                    r_en_un4.next = 1
#                    if i_w_muxed_wr_addr == i_w_muxed_wr_addr.max - 1:
#                        r_en_un4.next = 0
#                        r_state.next = State.IDLE

        fsm_union = self.FSM_UNION4(i_clk, r_en_un4, o_w_union4_rd_addr, o_w_dualram_data, o_w_union4_we, o_w_union4_data,o_w_union4_wr_addr, FFT_SIZE, 4, r_log_val)

        @always_comb
        def f_o_active_write():
            if r_log_val == 8:
                o_active_write.next = i_w_muxed_we
            else:
                o_active_write.next = 0
        @always_comb
        def f_o_data():
            o_dout.next = o_w_union4_data
        @always_comb
        def f_o_addr():
            o_addr.next = o_w_union4_wr_addr
        
        inst = instances()
        return inst#f_o_active_write, f_o_data,f_o_addr, fsm_input, f_muxed_data, f_muxed_rd_addr, f_muxed_wr_addr, f_muxed_we, seq, fsm_dualram, fsm_union#, fsm_bf2
    def tb(self):
        i_clk = Signal(bool(0))
        i_rst = Signal(bool(0))
        i_en = Signal(bool(0))
        o_union = Signal(bool(0))
        FFT_SIZE = 256
        o_addr = Signal(intbv(0,0,FFT_SIZE))
        #i_write_adr = Signal(modbv(0,0,2**N))
        
        cnt = Signal(intbv(0,0,512))
        i_new_data = Signal(bool(0))
        M = 16 - 1
        M_o = 32-1
        i_data = Signal(intbv(0, -2**M, 2**M))
        o_data = Signal(intbv(0, -2**M_o, 2**M_o))
        test_data = tuple([int(np.around(8192.0*np.cos(2*np.pi*1*i/FFT_SIZE))) for i in range(FFT_SIZE)]) ##TEST
#        print [hex(test_data[bit_reverse_addr[i]]) for i in range(FFT_SIZE)]

        #test_data = tuple([int(16384.0*np.cos(2*np.pi*i/FFT_SIZE)) << 16 for i in range(FFT_SIZE)]) ##TEST
#        uut = self.Shift_reg(i_clk, i_rst, i_en, i_new_data, i_data, o_union, o_data, o_addr, 16)
        FFT_in = self.FFT_FSM(i_clk, i_rst, i_en, i_new_data, i_data, o_data, o_union, o_addr,FFT_SIZE)

        i = Signal(modbv(0,0,FFT_SIZE))
        @instance
        def clk_gen():
            while True:
                yield delay(1)
                i_clk.next = not i_clk

        @instance
        def control():
            yield i_clk.posedge
            i_rst.next = 1
            yield i_clk.posedge
            yield i_clk.posedge
            yield i_clk.posedge
            i_data.next = 0
            i_en.next = 1
            i_rst.next = 0
            while True:
                yield i_clk.posedge
                cnt.next = (cnt + 1 ) % 512
                if cnt.next == 0:
                    i_new_data.next = 1
                    i_data.next = test_data[i]
                    i.next = i + 1
                else:
                    i_new_data.next = 0

        return clk_gen, control, FFT_in#, FFT_in

        
my_class = my_fft()
#
#
#tb = my_class.tb
#inst = traceSignals(tb)
#sim = Simulation(inst)
#sim.run(150000)

i_clk = Signal(bool(0))
i_rst = Signal(bool(0))
i_en = Signal(bool(0))
i_new_data = Signal(bool(0))
o_union = Signal(bool(0))
n = 16
n_signed = n - 1
N = 32
N_signed = N - 1
i_data = Signal(intbv(0, -2**n_signed, 2**n_signed))
o_data = Signal(intbv(0, -2**N_signed, 2**N_signed))
FFT_size = 256
o_addr = Signal(intbv(0,0,FFT_size))
toVerilog(my_class.FFT_FSM, i_clk, i_rst, i_en, i_new_data, i_data, o_data, o_union, o_addr,FFT_size)
#toVerilog(my_class.tb)
