from __future__ import annotations

import streamlit as st


def apply_theme() -> None:
    st.markdown(
        """
        <style>
        :root {
          --ms-bg: #08131f;
          --ms-bg-elevated: #0e2133;
          --ms-bg-panel: #132a3d;
          --ms-bg-input: #102638;
          --ms-bg-hover: #18364b;
          --ms-border: #28455b;
          --ms-border-strong: #37647d;
          --ms-text: #f5f8fa;
          --ms-text-secondary: #b2c2cf;
          --ms-text-muted: #8297a8;
          --ms-text-on-accent: #04151c;
          --ms-primary: #20c7c9;
          --ms-primary-hover: #43d7d3;
          --ms-success: #64e3a1;
          --ms-warning: #f3b45b;
          --ms-error: #ef6a72;
          --ms-focus: rgba(32, 199, 201, 0.35);
          --ms-accent: #ff4fa7;
          --ms-radius-sm: 6px;
          --ms-radius-md: 10px;
          --ms-radius-lg: 14px;

          --ms-ink: var(--ms-text);
          --ms-muted: var(--ms-text-secondary);
          --ms-line: var(--ms-border);
          --ms-panel: rgba(14, 33, 51, 0.82);
          --ms-teal: var(--ms-primary);
          --ms-mint: #64E3C4;
          --ms-deep: #06111d;
          --ms-green: var(--ms-success);
          --ms-amber: var(--ms-warning);
        }
        *,
        *::before,
        *::after {
          box-sizing: border-box;
        }
        .stApp {
          background:
            radial-gradient(circle at 20% 10%, rgba(32, 199, 201, 0.15), transparent 28rem),
            radial-gradient(circle at 85% 5%, rgba(255, 79, 167, 0.14), transparent 25rem),
            linear-gradient(140deg, var(--ms-bg) 0%, #0b1f2d 54%, #102b2e 100%);
          color: var(--ms-text);
        }
        html,
        body,
        #root,
        .stApp,
        [data-testid="stAppViewContainer"] {
          min-height: 100vh !important;
          overflow-y: auto !important;
        }
        .stApp > header,
        [data-testid="stHeader"],
        header[data-testid="stHeader"] {
          display: block !important;
          height: 3rem !important;
          min-height: 3rem !important;
          opacity: 1 !important;
          pointer-events: auto !important;
          background: transparent !important;
          background-color: transparent !important;
          box-shadow: none !important;
          overflow: visible !important;
          z-index: 999998 !important;
        }
        [data-testid="stHeader"]::before,
        [data-testid="stHeader"]::after {
          display: none !important;
          background: transparent !important;
        }
        [data-testid="stDecoration"],
        [data-testid="stStatusWidget"] {
          display: none !important;
        }
        [data-testid="stToolbar"] {
          background: transparent !important;
          pointer-events: auto !important;
        }
        [data-testid="collapsedControl"],
        [data-testid="stSidebarCollapsedControl"],
        button[data-testid="stBaseButton-headerNoPadding"],
        button[title="Open sidebar"],
        button[aria-label="Open sidebar"],
        button[title="Close sidebar"],
        button[aria-label="Close sidebar"] {
          display: flex !important;
          visibility: visible !important;
          opacity: 1 !important;
          pointer-events: auto !important;
          background: rgba(14, 33, 51, 0.92) !important;
          border: 1px solid rgba(40, 213, 245, 0.45) !important;
          border-radius: 8px !important;
          color: var(--ms-ink) !important;
          z-index: 999999 !important;
        }
        [data-testid="collapsedControl"] *,
        [data-testid="stSidebarCollapsedControl"] *,
        button[data-testid="stBaseButton-headerNoPadding"] *,
        button[title="Open sidebar"] *,
        button[aria-label="Open sidebar"] * {
          color: #ffffff !important;
          fill: #ffffff !important;
          stroke: #ffffff !important;
        }
        [data-testid="collapsedControl"] svg,
        [data-testid="stSidebarCollapsedControl"] svg,
        button[data-testid="stBaseButton-headerNoPadding"] svg,
        button[title="Open sidebar"] svg,
        button[aria-label="Open sidebar"] svg,
        button[title="Close sidebar"] svg,
        button[aria-label="Close sidebar"] svg {
          color: #ffffff !important;
          fill: #ffffff !important;
          stroke: #ffffff !important;
        }
        [data-testid="stSidebar"] {
          background: rgba(5, 16, 28, 0.98);
          border-right: 1px solid var(--ms-line);
          scrollbar-color: rgba(40, 213, 245, 0.78) rgba(8, 19, 31, 0.65);
          scrollbar-width: thin;
        }
        [data-testid="stSidebarContent"] {
          overflow-y: auto !important;
          max-height: 100vh !important;
          padding-top: 1.25rem;
          scrollbar-color: rgba(40, 213, 245, 0.78) rgba(8, 19, 31, 0.65);
          scrollbar-width: thin;
        }
        [data-testid="stSidebar"]::-webkit-scrollbar,
        [data-testid="stSidebarContent"]::-webkit-scrollbar,
        [data-testid="stAppViewContainer"]::-webkit-scrollbar,
        body::-webkit-scrollbar {
          width: 10px;
          height: 10px;
        }
        [data-testid="stSidebar"]::-webkit-scrollbar-track,
        [data-testid="stSidebarContent"]::-webkit-scrollbar-track,
        [data-testid="stAppViewContainer"]::-webkit-scrollbar-track,
        body::-webkit-scrollbar-track {
          background: rgba(8, 19, 31, 0.72);
        }
        [data-testid="stSidebar"]::-webkit-scrollbar-thumb,
        [data-testid="stSidebarContent"]::-webkit-scrollbar-thumb,
        [data-testid="stAppViewContainer"]::-webkit-scrollbar-thumb,
        body::-webkit-scrollbar-thumb {
          background: linear-gradient(180deg, var(--ms-teal), var(--ms-accent));
          border: 2px solid rgba(8, 19, 31, 0.72);
          border-radius: 999px;
        }
        [data-testid="stSidebar"]::-webkit-scrollbar-thumb:hover,
        [data-testid="stSidebarContent"]::-webkit-scrollbar-thumb:hover,
        [data-testid="stAppViewContainer"]::-webkit-scrollbar-thumb:hover,
        body::-webkit-scrollbar-thumb:hover {
          background: var(--ms-mint, #64E3C4);
        }
        [data-testid="stSidebar"] {
          color: var(--ms-text);
        }
        [data-testid="stSidebar"] p,
        [data-testid="stSidebar"] label,
        [data-testid="stSidebar"] span {
          color: inherit;
        }
        [data-testid="stSidebar"] img {
          background: rgba(255,255,255,0.92);
          border: 1px solid rgba(255,255,255,0.18);
          border-radius: 8px;
          padding: 0.8rem;
          max-height: 21rem;
          object-fit: contain;
        }
        .block-container {
          padding-top: 3.4rem;
          padding-bottom: 2rem;
          max-width: 1500px;
        }
        h1, h2, h3 {
          letter-spacing: 0;
        }
        h2 {
          color: var(--ms-ink);
        }
        div[data-testid="stMetric"] {
          background: var(--ms-bg-elevated);
          border: 1px solid var(--ms-line);
          border-radius: 8px;
          padding: 0.75rem 0.9rem;
        }
        div[data-testid="stMetric"] label,
        div[data-testid="stMetric"] [data-testid="stMetricValue"] {
          color: var(--ms-text);
        }
        div[data-testid="stMetric"] label {
          color: var(--ms-text-secondary) !important;
        }
        .ms-header {
          border-bottom: 1px solid var(--ms-line);
          padding: 0.4rem 0 1.35rem 0;
          margin-bottom: 1.1rem;
          display: grid;
          grid-template-columns: minmax(0, 1fr) 7.5rem;
          gap: 1.25rem;
          align-items: center;
        }
        .ms-kicker {
          color: var(--ms-teal);
          font-size: 0.82rem;
          font-weight: 700;
          text-transform: uppercase;
        }
        .ms-brand {
          color: var(--ms-accent);
          font-size: clamp(3.6rem, 8vw, 7.4rem);
          font-weight: 900;
          line-height: 0.88;
          margin: 0.18rem 0 0.45rem;
        }
        .ms-soft {
          color: var(--ms-muted);
        }
        .ms-sidebar-title {
          color: var(--ms-accent);
          font-size: 1.55rem;
          font-weight: 900;
          line-height: 1;
          margin: -0.4rem 0 0.8rem;
        }
        .ms-topbar {
          display: grid;
          grid-template-columns: minmax(18rem, 1fr) auto auto;
          gap: 1rem;
          align-items: center;
          border: 1px solid var(--ms-line);
          background: rgba(5, 16, 28, 0.72);
          border-radius: 8px;
          padding: 0.72rem 1rem;
          margin-bottom: 1rem;
          box-shadow: 0 18px 48px rgba(0,0,0,0.22);
        }
        .ms-top-brand {
          display: flex;
          gap: 0.85rem;
          align-items: center;
          min-width: 0;
        }
        .ms-top-brand img {
          width: 3.8rem;
          height: 3.8rem;
          object-fit: contain;
          background: rgba(255,255,255,0.93);
          border-radius: 8px;
          padding: 0.28rem;
        }
        .ms-top-title {
          color: var(--ms-accent);
          font-size: clamp(1.9rem, 3.2vw, 3.35rem);
          font-weight: 950;
          line-height: 0.92;
        }
        .ms-top-subtitle {
          color: var(--ms-muted);
          font-size: 0.88rem;
          margin-top: 0.2rem;
          white-space: nowrap;
          overflow: hidden;
          text-overflow: ellipsis;
          max-width: 34rem;
        }
        .ms-top-pills {
          display: flex;
          gap: 0.5rem;
          flex-wrap: wrap;
          justify-content: flex-end;
        }
        .ms-top-pills span,
        .ms-ready {
          border: 1px solid var(--ms-line);
          border-radius: 8px;
          background: rgba(255,255,255,0.06);
          padding: 0.54rem 0.78rem;
          color: var(--ms-ink);
          font-weight: 700;
        }
        .ms-top-pills span:nth-child(2) {
          color: var(--ms-teal);
          border-color: rgba(40, 213, 245, 0.4);
        }
        .ms-top-pills span:nth-child(3) {
          color: var(--ms-green);
          border-color: rgba(91, 234, 137, 0.35);
        }
        .ms-ready {
          display: flex;
          gap: 0.45rem;
          align-items: center;
          justify-content: center;
        }
        .ms-ready span {
          width: 0.62rem;
          height: 0.62rem;
          border-radius: 50%;
          background: var(--ms-green);
          box-shadow: 0 0 18px rgba(91, 234, 137, 0.8);
        }
        .ms-steps {
          margin: 1.4rem 0 1rem;
        }
        .ms-step {
          display: flex;
          align-items: center;
          gap: 0.7rem;
          padding: 0.58rem 0.2rem;
          color: var(--ms-muted);
        }
        .ms-step span {
          width: 1.82rem;
          height: 1.82rem;
          border-radius: 50%;
          border: 1px solid rgba(40, 213, 245, 0.45);
          background: rgba(255,255,255,0.04);
        }
        .ms-step p {
          margin: 0;
          font-weight: 700;
        }
        .ms-step.active {
          color: var(--ms-deep);
          background: linear-gradient(135deg, var(--ms-teal), #7ce7f4);
          border-radius: 8px;
          padding-left: 0.7rem;
          margin: 0.25rem 0;
        }
        .ms-step.active span {
          border-color: var(--ms-deep);
        }
        .ms-step.done span {
          background: var(--ms-green);
          border-color: var(--ms-green);
        }
        .ms-control-title,
        .ms-panel-title {
          color: var(--ms-ink);
          font-size: 1.22rem;
          font-weight: 850;
          margin: 0.25rem 0 0.75rem;
        }
        div[data-testid="stTabs"] {
          border: 1px solid var(--ms-line);
          border-radius: 8px;
          background: rgba(14, 33, 51, 0.66);
          padding: 0.75rem;
        }
        div[data-testid="stTabs"] button {
          color: var(--ms-muted) !important;
        }
        div[data-testid="stTabs"] button p {
          color: inherit !important;
          font-weight: 750;
        }
        div[data-testid="stTabs"] button[aria-selected="true"] {
          color: var(--ms-deep) !important;
          background: linear-gradient(135deg, var(--ms-primary), var(--ms-primary-hover));
          border-radius: 8px;
        }
        div[data-testid="stRadio"] label,
        div[data-testid="stCheckbox"] label,
        div[data-testid="stToggle"] label {
          color: var(--ms-text) !important;
        }
        div[data-testid="stRadio"] label:hover,
        div[data-testid="stCheckbox"] label:hover,
        div[data-testid="stToggle"] label:hover {
          color: var(--ms-primary-hover) !important;
        }
        div[data-testid="stProgress"] > div > div {
          background-color: var(--ms-primary) !important;
        }
        [data-testid="stWidgetLabel"] p,
        [data-testid="stWidgetLabel"] label,
        [data-testid="stWidgetLabel"] span {
          color: #dce8ef !important;
          font-weight: 700;
        }
        div[data-testid="stTextInput"] input,
        div[data-testid="stTextArea"] textarea,
        div[data-testid="stNumberInput"] input {
          background-color: var(--ms-bg-input) !important;
          color: var(--ms-text) !important;
          border: 1px solid var(--ms-border) !important;
          border-radius: var(--ms-radius-sm) !important;
          box-shadow: none !important;
        }
        div[data-testid="stTextInput"] input:hover,
        div[data-testid="stTextArea"] textarea:hover,
        div[data-testid="stNumberInput"] input:hover {
          border-color: var(--ms-border-strong) !important;
          background-color: var(--ms-bg-hover) !important;
        }
        div[data-testid="stTextInput"] input:focus,
        div[data-testid="stTextArea"] textarea:focus,
        div[data-testid="stNumberInput"] input:focus {
          border-color: var(--ms-primary) !important;
          box-shadow: 0 0 0 0.18rem var(--ms-focus) !important;
        }
        div[data-baseweb="select"] > div,
        div[data-testid="stMultiSelect"] div[data-baseweb="select"] > div {
          background-color: var(--ms-bg-input) !important;
          color: var(--ms-text) !important;
          border-color: var(--ms-border) !important;
          border-radius: var(--ms-radius-sm) !important;
        }
        div[data-baseweb="select"] span,
        div[data-baseweb="select"] svg {
          color: var(--ms-text-secondary) !important;
          fill: var(--ms-text-secondary) !important;
        }
        div[data-baseweb="popover"],
        div[data-baseweb="menu"],
        ul[role="listbox"] {
          background: var(--ms-bg-elevated) !important;
          color: var(--ms-text) !important;
          border: 1px solid var(--ms-border-strong) !important;
          border-radius: var(--ms-radius-md) !important;
          box-shadow: 0 18px 50px rgba(0,0,0,0.42) !important;
        }
        li[role="option"],
        div[role="option"] {
          background: transparent !important;
          color: var(--ms-text) !important;
        }
        li[role="option"]:hover,
        div[role="option"]:hover,
        li[role="option"][aria-selected="true"],
        div[role="option"][aria-selected="true"] {
          background: var(--ms-bg-hover) !important;
          color: var(--ms-text) !important;
        }
        div[data-baseweb="tag"] {
          background: rgba(32, 199, 201, 0.18) !important;
          border: 1px solid rgba(32, 199, 201, 0.38) !important;
          color: var(--ms-text) !important;
        }
        div[data-baseweb="tag"] span,
        div[data-baseweb="tag"] svg {
          color: var(--ms-text) !important;
          fill: var(--ms-text) !important;
        }
        input::placeholder,
        textarea::placeholder {
          color: var(--ms-text-muted) !important;
        }
        [data-testid="stFileUploader"] {
          margin-bottom: 0.8rem;
        }
        [data-testid="stFileUploader"] section,
        [data-testid="stFileUploaderDropzone"] {
          background: var(--ms-bg-input) !important;
          border: 1px dashed rgba(40, 213, 245, 0.36) !important;
          border-radius: 8px !important;
          min-height: 4.4rem;
        }
        [data-testid="stFileUploader"] section p,
        [data-testid="stFileUploader"] section span,
        [data-testid="stFileUploader"] section small {
          color: var(--ms-text-secondary) !important;
        }
        [data-testid="stFileUploader"] button {
          background: rgba(40, 213, 245, 0.14) !important;
          border: 1px solid rgba(40, 213, 245, 0.42) !important;
          color: var(--ms-ink) !important;
        }
        [data-testid="stExpander"] {
          border: 1px solid var(--ms-line);
          background: rgba(14,33,51,0.62);
          border-radius: 8px;
        }
        [data-testid="stExpander"] summary,
        [data-testid="stExpander"] summary p {
          color: #dce8ef !important;
          font-weight: 750;
        }
        [data-testid="stMarkdownContainer"] code,
        pre {
          color: var(--ms-text) !important;
          background: #07111c !important;
          border-color: var(--ms-border) !important;
        }
        .ms-canvas {
          position: relative;
          min-height: 520px;
          border: 1px solid var(--ms-line);
          border-radius: 8px;
          overflow: hidden;
          background:
            radial-gradient(circle at 44% 28%, rgba(40, 213, 245, 0.22), transparent 18rem),
            radial-gradient(circle at 70% 30%, rgba(255, 79, 167, 0.16), transparent 18rem),
            linear-gradient(180deg, rgba(5,16,28,0.78), rgba(2,8,15,0.92));
          box-shadow: inset 0 0 90px rgba(0,0,0,0.38), 0 18px 50px rgba(0,0,0,0.28);
        }
        .ms-legend {
          position: absolute;
          top: 1.1rem;
          left: 1.15rem;
          display: grid;
          gap: 0.45rem;
          z-index: 3;
        }
        .ms-legend span {
          display: flex;
          align-items: center;
          gap: 0.5rem;
          color: var(--ms-ink);
          font-weight: 700;
        }
        .ms-legend b {
          display: inline-block;
          width: 0.75rem;
          height: 0.75rem;
          border-radius: 50%;
        }
        .ms-legend .cyan { background: var(--ms-teal); }
        .ms-legend .amber { background: var(--ms-amber); }
        .ms-legend .silver { background: #aeb8c1; }
        .ms-protein-shape {
          position: absolute;
          left: 21%;
          right: 10%;
          top: 8%;
          height: 50%;
          filter: drop-shadow(0 28px 38px rgba(0,0,0,0.38));
        }
        .ms-protein-shape i {
          position: absolute;
          display: block;
          border-radius: 50%;
          background: radial-gradient(circle at 28% 25%, #a9f6ff, var(--ms-teal) 45%, #176b87);
          opacity: 0.94;
        }
        .ms-protein-shape i:nth-child(1) { width: 16rem; height: 13rem; left: 0; top: 2.2rem; }
        .ms-protein-shape i:nth-child(2) { width: 13rem; height: 12rem; left: 11rem; top: 0.5rem; background: radial-gradient(circle at 30% 25%, #c9bcff, #6967d8 50%, #343681); }
        .ms-protein-shape i:nth-child(3) { width: 12rem; height: 11rem; left: 5.5rem; top: 11rem; background: radial-gradient(circle at 30% 25%, #8bd8ff, #357de1 55%, #174b92); }
        .ms-protein-shape i:nth-child(4) { width: 13rem; height: 10.5rem; left: 21rem; top: 11.4rem; background: radial-gradient(circle at 30% 25%, #c5fff5, #61c8c0 52%, #286c69); }
        .ms-protein-shape i:nth-child(5) { width: 8rem; height: 8rem; left: 27rem; top: 4rem; background: radial-gradient(circle at 30% 25%, #ff9fca, var(--ms-accent) 55%, #9b1f61); }
        .ms-protein-shape i:nth-child(6) { width: 7rem; height: 7rem; left: 14.5rem; top: 8.5rem; background: radial-gradient(circle at 30% 25%, #ffffff, #76e4ef 52%, #2d7e95); }
        .ms-surface-grid {
          position: absolute;
          left: -4%;
          right: -4%;
          bottom: -7%;
          height: 34%;
          background:
            radial-gradient(circle, rgba(210, 222, 231, 0.36) 0 0.48rem, transparent 0.52rem) 0 0 / 2.1rem 1.35rem,
            radial-gradient(circle, rgba(130, 148, 160, 0.38) 0 0.46rem, transparent 0.52rem) 1.05rem 0.68rem / 2.1rem 1.35rem,
            linear-gradient(180deg, rgba(255,255,255,0.12), rgba(0,0,0,0.24));
          transform: perspective(580px) rotateX(58deg);
          transform-origin: bottom;
        }
        .ms-anchor-callout {
          position: absolute;
          bottom: 34%;
          border: 1px solid rgba(255, 185, 76, 0.55);
          background: rgba(5, 16, 28, 0.82);
          color: var(--ms-amber);
          border-radius: 8px;
          padding: 0.45rem 0.65rem;
          font-weight: 800;
          z-index: 4;
        }
        .ms-anchor-callout.left { left: 9%; }
        .ms-anchor-callout.right { right: 8%; color: var(--ms-teal); border-color: rgba(40,213,245,0.45); }
        .ms-axis {
          position: absolute;
          left: 2rem;
          bottom: 2rem;
          color: var(--ms-ink);
          font-weight: 900;
        }
        .ms-spec-grid {
          display: grid;
          grid-template-columns: repeat(4, minmax(0, 1fr));
          gap: 0.8rem;
          margin: 1rem 0;
        }
        .ms-spec {
          border: 1px solid var(--ms-line);
          border-radius: 8px;
          background: rgba(10, 26, 42, 0.72);
          padding: 0.95rem;
        }
        .ms-spec span {
          display: block;
          color: var(--ms-muted);
          font-size: 0.82rem;
          margin-bottom: 0.32rem;
        }
        .ms-spec strong {
          color: var(--ms-ink);
          font-size: 1rem;
        }
        .ms-command-bar {
          border: 1px solid var(--ms-line);
          border-radius: 8px;
          background: rgba(5, 16, 28, 0.72);
          padding: 0.6rem;
        }
        .ms-empty-preview {
          min-height: 430px;
          border: 1px dashed rgba(40,213,245,0.32);
          border-radius: 8px;
          background:
            radial-gradient(circle at center, rgba(40,213,245,0.12), transparent 18rem),
            rgba(5,16,28,0.52);
          display: grid;
          place-content: center;
          text-align: center;
          padding: 2rem;
        }
        .ms-empty-preview strong {
          color: var(--ms-ink);
          font-size: 1.4rem;
        }
        .ms-empty-preview span {
          color: var(--ms-muted);
          max-width: 32rem;
          margin-top: 0.45rem;
        }
        .ms-surface-preview {
          margin-top: 1rem;
          border: 1px solid var(--ms-line);
          border-radius: 8px;
          background: rgba(5,16,28,0.54);
          padding: 1rem;
        }
        .ms-surface-preview span {
          display: block;
          color: var(--ms-muted);
          font-weight: 700;
          margin-bottom: 0.8rem;
        }
        .ms-surface-preview div {
          height: 9rem;
          border-radius: 8px;
          background:
            radial-gradient(circle, rgba(40,213,245,0.72) 0 0.32rem, transparent 0.34rem) 0 0 / 1.35rem 1.05rem,
            radial-gradient(circle, rgba(100,227,196,0.62) 0 0.30rem, transparent 0.34rem) 0.68rem 0.52rem / 1.35rem 1.05rem,
            linear-gradient(180deg, rgba(255,255,255,0.06), rgba(0,0,0,0.18));
          transform: perspective(460px) rotateX(54deg);
          transform-origin: bottom;
        }
        .ms-quality {
          margin-top: 1rem;
          border: 1px solid var(--ms-line);
          border-radius: 8px;
          background: rgba(255,255,255,0.045);
          padding: 0.9rem;
        }
        .ms-quality strong {
          color: var(--ms-ink);
          display: block;
          margin-bottom: 0.6rem;
        }
        .ms-quality-row {
          display: flex;
          align-items: center;
          gap: 0.5rem;
          color: var(--ms-muted);
          padding: 0.25rem 0;
          font-weight: 700;
        }
        .ms-quality-row span {
          width: 0.68rem;
          height: 0.68rem;
          border-radius: 50%;
          background: var(--ms-muted);
        }
        .ms-quality-row.ok {
          color: var(--ms-green);
        }
        .ms-quality-row.ok span {
          background: var(--ms-green);
          box-shadow: 0 0 12px rgba(91,234,137,0.55);
        }
        .ms-hero-logo img {
          max-height: 12rem;
          object-fit: contain;
        }
        .ms-chip-row {
          display: flex;
          gap: 0.45rem;
          flex-wrap: wrap;
          margin-top: 0.8rem;
        }
        .ms-chip {
          border: 1px solid var(--ms-line);
          border-radius: 999px;
          color: var(--ms-deep);
          background: rgba(255,255,255,0.72);
          padding: 0.28rem 0.68rem;
          font-size: 0.82rem;
          font-weight: 650;
        }
        .ms-section {
          border: 1px solid var(--ms-line);
          border-radius: 8px;
          background: var(--ms-bg-elevated);
          padding: 1rem;
          margin-bottom: 1rem;
        }
        .ms-log-row {
          border-left: 3px solid var(--ms-accent);
          padding: 0.35rem 0 0.35rem 0.75rem;
          margin: 0.35rem 0;
          background: rgba(255,255,255,0.06);
          border-radius: 0 8px 8px 0;
        }
        .ms-log-label {
          color: var(--ms-teal);
          font-size: 0.78rem;
          font-weight: 750;
          text-transform: uppercase;
        }
        .ms-log-value {
          color: var(--ms-ink);
          font-size: 0.92rem;
          margin-top: 0.08rem;
        }
        .stButton > button, .stDownloadButton > button {
          border-radius: 8px;
          border: 1px solid var(--ms-border-strong);
          background: var(--ms-bg-elevated);
          color: var(--ms-text);
        }
        .stButton > button:hover,
        .stDownloadButton > button:hover {
          border-color: var(--ms-primary);
          background: var(--ms-bg-hover);
          color: var(--ms-text);
        }
        .stButton > button:focus,
        .stDownloadButton > button:focus,
        .stLinkButton > a:focus {
          box-shadow: 0 0 0 0.18rem var(--ms-focus) !important;
          outline: none !important;
        }
        .stButton > button:disabled,
        .stDownloadButton > button:disabled {
          border-color: rgba(40, 69, 91, 0.55) !important;
          background: rgba(16, 38, 56, 0.46) !important;
          color: var(--ms-text-muted) !important;
          cursor: not-allowed !important;
        }
        .stLinkButton > a,
        div[data-testid="stPopover"] button,
        .stDownloadButton > button {
          border-radius: 8px !important;
          border: 1px solid var(--ms-border-strong) !important;
          background: var(--ms-bg-elevated) !important;
          color: var(--ms-text) !important;
          font-weight: 800 !important;
        }
        .stLinkButton > a p,
        div[data-testid="stPopover"] button p,
        .stDownloadButton > button p {
          color: var(--ms-text) !important;
          overflow-wrap: normal !important;
          word-break: normal !important;
          white-space: normal !important;
        }
        .stButton > button[kind="primary"] {
          background: linear-gradient(135deg, var(--ms-primary), var(--ms-primary-hover));
          color: var(--ms-text-on-accent);
          border-color: var(--ms-primary);
          font-weight: 850;
        }
        div[data-testid="stAlert"] {
          border-radius: var(--ms-radius-md);
          border: 1px solid var(--ms-border);
          background: var(--ms-bg-elevated);
          color: var(--ms-text);
        }
        div[data-testid="stAlert"] p,
        div[data-testid="stAlert"] li {
          color: var(--ms-text) !important;
        }
        div[data-testid="stDataFrame"],
        div[data-testid="stDataEditor"] {
          border: 1px solid var(--ms-border);
          border-radius: var(--ms-radius-md);
          overflow: hidden;
        }
        div[data-testid="stDataFrame"] *,
        div[data-testid="stDataEditor"] * {
          color: inherit;
        }
        div[data-testid="stStatusWidget"],
        div[data-testid="stStatus"] {
          color: var(--ms-text) !important;
        }
        @media (max-width: 900px) {
          .ms-header {
            grid-template-columns: 1fr;
          }
          .ms-topbar {
            grid-template-columns: 1fr;
          }
          .ms-hero-logo {
            display: none;
          }
          .ms-spec-grid {
            grid-template-columns: 1fr 1fr;
          }
          .ms-canvas {
            min-height: 420px;
          }
        }
        @media (max-width: 1100px) {
          .block-container {
            padding-left: 1.25rem;
            padding-right: 1.25rem;
          }
          section[data-testid="stMain"] div[data-testid="stHorizontalBlock"] {
            display: block !important;
          }
          section[data-testid="stMain"] div[data-testid="stHorizontalBlock"] > div[data-testid="stColumn"] {
            flex: 1 1 100% !important;
            min-width: 100% !important;
            width: 100% !important;
            margin-bottom: 0.75rem;
          }
          .stLinkButton > a,
          div[data-testid="stPopover"] button,
          .stDownloadButton > button,
          .stButton > button {
            min-height: 2.75rem;
          }
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
