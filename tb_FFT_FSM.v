`timescale 1ns/10ps
module tb_FFT_FSM;

reg i_clk;
reg i_rst;
reg i_en;
reg i_new_data;
reg [15:0] i_data;
wire [31:0] o_dout;
wire o_active_write;
wire [7:0] o_addr;


FFT_FSM dut(
    i_clk,
    i_rst,
    i_en,
    i_new_data,
    i_data,
    o_dout,
    o_active_write,
    o_addr
);
initial begin
    i_clk = 0;
    i_rst = 1;
    i_en = 0;
    i_new_data = 0;
    i_data = 0;
    #100
    i_rst <= #20 0;
    i_en <= #20 1; 
end
always #10 i_clk <=  !i_clk;
endmodule
