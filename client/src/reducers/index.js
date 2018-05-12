import { combineReducers } from "redux";
import context from "./context";
import jobs from "./jobs";

export default combineReducers({ context, jobs });
