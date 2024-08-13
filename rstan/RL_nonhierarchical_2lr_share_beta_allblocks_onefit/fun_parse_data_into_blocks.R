#function to parse vol_train data for R

fun_parse_data_into_blocks<- function(datalist){
  datain<-datalist$data
  subnum<-datalist$subnum
  blktype<-datalist$blktype
  datain$button_press<-abs(datain$button_press-2)#change 1(shape1chosen) to1 (shape1chosen);2(shape2chosen) to 0(shape1 not chosen)

  block1_trials<-datain[which(datain$trial_number<81),]
  block2_trials<-datain[which(datain$trial_number>80 & datain$trial_number<161),]
  block3_trials<-datain[which(datain$trial_number>160),]
  
  
  #info is 1 if outcome associated with shape 1, 0 shape 2
  #choices are 1 for chosing shape 1, 0 shape 2
  block1_wins<-block1_trials$shape1_win
  block1_loss<-block1_trials$shape1_loss
  block1_choice<-block1_trials$button_press
  block1_RT<-block1_trials$reaction_time
  block2_wins<-block2_trials$shape1_win
  block2_loss<-block2_trials$shape1_loss
  block2_choice<-block2_trials$button_press
  block2_RT<-block2_trials$reaction_time
  block3_wins<-block3_trials$shape1_win
  block3_loss<-block3_trials$shape1_loss
  block3_choice<-block3_trials$button_press
  block3_RT<-block3_trials$reaction_time
  
  return(ret_list<-list("data" = data.frame(block1_wins, block1_loss, block1_choice,block1_RT,block2_wins, block2_loss, block2_choice, block2_RT,
                                            block3_wins, block3_loss, block3_choice,block3_RT), 
                                            "subnum"=subnum,"blktype"=blktype))
}