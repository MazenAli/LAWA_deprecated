namespace lawa {

template <typename RB_Model, typename TruthModel>
RB_Base<RB_Model,TruthModel>::RB_Base(RB_Model& _rb_system, TruthModel& _rb_truth)
 : rb_system(_rb_system), rb_truth(_rb_truth){}

} // namespace lawa
