load disney;
addpath('../helpers');

Draw_Bollinger_Bands( 400:600, 20, 1, dis_HIGH, dis_LOW, dis_CLOSE, dis_OPEN, false );
Draw_Bollinger_Bands( 400:600, 20, 2, dis_HIGH, dis_LOW, dis_CLOSE, dis_OPEN, false );
Draw_Bollinger_Bands( 400:600, 40, 1.5, dis_HIGH, dis_LOW, dis_CLOSE, dis_OPEN, false );
Draw_Bollinger_Bands( 400:600, 40, 2.5, dis_HIGH, dis_LOW, dis_CLOSE, dis_OPEN, false );