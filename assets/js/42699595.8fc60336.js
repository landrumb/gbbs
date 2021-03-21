(window.webpackJsonp=window.webpackJsonp||[]).push([[10],{120:function(e,a,t){"use strict";t.d(a,"a",(function(){return i})),t.d(a,"b",(function(){return h}));var n=t(0),s=t.n(n);function m(e,a,t){return a in e?Object.defineProperty(e,a,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[a]=t,e}function p(e,a){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var n=Object.getOwnPropertySymbols(e);a&&(n=n.filter((function(a){return Object.getOwnPropertyDescriptor(e,a).enumerable}))),t.push.apply(t,n)}return t}function r(e){for(var a=1;a<arguments.length;a++){var t=null!=arguments[a]?arguments[a]:{};a%2?p(Object(t),!0).forEach((function(a){m(e,a,t[a])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):p(Object(t)).forEach((function(a){Object.defineProperty(e,a,Object.getOwnPropertyDescriptor(t,a))}))}return e}function c(e,a){if(null==e)return{};var t,n,s=function(e,a){if(null==e)return{};var t,n,s={},m=Object.keys(e);for(n=0;n<m.length;n++)t=m[n],a.indexOf(t)>=0||(s[t]=e[t]);return s}(e,a);if(Object.getOwnPropertySymbols){var m=Object.getOwnPropertySymbols(e);for(n=0;n<m.length;n++)t=m[n],a.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(s[t]=e[t])}return s}var b=s.a.createContext({}),l=function(e){var a=s.a.useContext(b),t=a;return e&&(t="function"==typeof e?e(a):r(r({},a),e)),t},i=function(e){var a=l(e.components);return s.a.createElement(b.Provider,{value:a},e.children)},N={inlineCode:"code",wrapper:function(e){var a=e.children;return s.a.createElement(s.a.Fragment,{},a)}},o=s.a.forwardRef((function(e,a){var t=e.components,n=e.mdxType,m=e.originalType,p=e.parentName,b=c(e,["components","mdxType","originalType","parentName"]),i=l(t),o=n,h=i["".concat(p,".").concat(o)]||i[o]||N[o]||m;return t?s.a.createElement(h,r(r({ref:a},b),{},{components:t})):s.a.createElement(h,r({ref:a},b))}));function h(e,a){var t=arguments,n=a&&a.mdxType;if("string"==typeof e||n){var m=t.length,p=new Array(m);p[0]=o;var r={};for(var c in a)hasOwnProperty.call(a,c)&&(r[c]=a[c]);r.originalType=e,r.mdxType="string"==typeof e?e:n,p[1]=r;for(var b=2;b<m;b++)p[b]=t[b];return s.a.createElement.apply(null,p)}return s.a.createElement.apply(null,t)}o.displayName="MDXCreateElement"},80:function(e,a,t){"use strict";t.r(a),t.d(a,"frontMatter",(function(){return p})),t.d(a,"metadata",(function(){return r})),t.d(a,"toc",(function(){return c})),t.d(a,"default",(function(){return l}));var n=t(3),s=t(7),m=(t(0),t(120)),p={id:"maximal_independent_set",title:"Maximal Independent Set"},r={unversionedId:"benchmarks/covering/maximal_independent_set",id:"benchmarks/covering/maximal_independent_set",isDocsHomePage:!1,title:"Maximal Independent Set",description:"Problem Specification",source:"@site/docs/benchmarks/covering/maximal_independent_set.md",slug:"/benchmarks/covering/maximal_independent_set",permalink:"/gbbs/docs/benchmarks/covering/maximal_independent_set",version:"current",sidebar:"docs",previous:{title:"Strongly Connected Components",permalink:"/gbbs/docs/benchmarks/connectivity/strongly_connected_components"},next:{title:"Maximal Matching",permalink:"/gbbs/docs/benchmarks/covering/maximal_matching"}},c=[{value:"Problem Specification",id:"problem-specification",children:[]},{value:"Algorithm Implementations",id:"algorithm-implementations",children:[]},{value:"Cost Bounds",id:"cost-bounds",children:[]},{value:"Compiling and Running",id:"compiling-and-running",children:[]},{value:"References",id:"references",children:[]}],b={toc:c};function l(e){var a=e.components,t=Object(s.a)(e,["components"]);return Object(m.b)("wrapper",Object(n.a)({},b,t,{components:a,mdxType:"MDXLayout"}),Object(m.b)("h2",{id:"problem-specification"},"Problem Specification"),Object(m.b)("h4",{id:"input"},"Input"),Object(m.b)("p",null,Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"G"),Object(m.b)("mo",{parentName:"mrow"},"="),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(m.b)("mi",{parentName:"mrow"},"V"),Object(m.b)("mo",{parentName:"mrow",separator:"true"},","),Object(m.b)("mi",{parentName:"mrow"},"E"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},")")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"G=(V, E)")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"G"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}}),Object(m.b)("span",{parentName:"span",className:"mrel"},"="),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}})),Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mopen"},"("),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.22222em"}},"V"),Object(m.b)("span",{parentName:"span",className:"mpunct"},","),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.16666666666666666em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.05764em"}},"E"),Object(m.b)("span",{parentName:"span",className:"mclose"},")"))))),", an undirected graph."),Object(m.b)("h4",{id:"output"},"Output"),Object(m.b)("p",null,Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"U"),Object(m.b)("mo",{parentName:"mrow"},"\u2286"),Object(m.b)("mi",{parentName:"mrow"},"V")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"U \\subseteq V")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.8193em",verticalAlign:"-0.13597em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.10903em"}},"U"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}}),Object(m.b)("span",{parentName:"span",className:"mrel"},"\u2286"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}})),Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.22222em"}},"V"))))),", a set of vertices such that no two vertices in ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"U")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"U")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.10903em"}},"U")))))," are neighbors and all vertices in ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"V"),Object(m.b)("mo",{parentName:"mrow"},"\u2216"),Object(m.b)("mi",{parentName:"mrow"},"U")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"V \\setminus U")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.22222em"}},"V"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2222222222222222em"}}),Object(m.b)("span",{parentName:"span",className:"mbin"},"\u2216"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2222222222222222em"}})),Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.10903em"}},"U")))))," have a neighbor in ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"U")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"U")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.10903em"}},"U"))))),"."),Object(m.b)("h2",{id:"algorithm-implementations"},"Algorithm Implementations"),Object(m.b)("p",null,"The maximal independent set problem is to compute a subset of vertices\n",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"U")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"U")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.10903em"}},"U")))))," such that no two vertices in ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"U")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"U")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.10903em"}},"U")))))," are neighbors, and all vertices\nin ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"V"),Object(m.b)("mo",{parentName:"mrow"},"\u2216"),Object(m.b)("mi",{parentName:"mrow"},"U")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"V \\setminus U")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.22222em"}},"V"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2222222222222222em"}}),Object(m.b)("span",{parentName:"span",className:"mbin"},"\u2216"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2222222222222222em"}})),Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.10903em"}},"U")))))," have a neighbor in ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"U")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"U")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.10903em"}},"U"))))),".  Maximal independent set\n(MIS) and maximal matching (MM) are easily solved in linear work\nsequentially using greedy algorithms."),Object(m.b)("p",null,"In GBBS, we implement the rootset-based algorithm for MIS from\nBlelloch et al.","[1]"," which runs in ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"O"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(m.b)("mi",{parentName:"mrow"},"m"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},")")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"O(m)")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.02778em"}},"O"),Object(m.b)("span",{parentName:"span",className:"mopen"},"("),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"m"),Object(m.b)("span",{parentName:"span",className:"mclose"},")")))))," work\nand ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"O"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(m.b)("msup",{parentName:"mrow"},Object(m.b)("mo",{parentName:"msup"},Object(m.b)("mi",{parentName:"mo"},"log"),Object(m.b)("mo",{parentName:"mo"},"\u2061")),Object(m.b)("mn",{parentName:"msup"},"2")),Object(m.b)("mi",{parentName:"mrow"},"n"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},")")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"O(\\log^2 n)")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1.148448em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.02778em"}},"O"),Object(m.b)("span",{parentName:"span",className:"mopen"},"("),Object(m.b)("span",{parentName:"span",className:"mop"},Object(m.b)("span",{parentName:"span",className:"mop"},"lo",Object(m.b)("span",{parentName:"span",style:{marginRight:"0.01389em"}},"g")),Object(m.b)("span",{parentName:"span",className:"msupsub"},Object(m.b)("span",{parentName:"span",className:"vlist-t"},Object(m.b)("span",{parentName:"span",className:"vlist-r"},Object(m.b)("span",{parentName:"span",className:"vlist",style:{height:"0.8984479999999999em"}},Object(m.b)("span",{parentName:"span",style:{top:"-3.1473400000000002em",marginRight:"0.05em"}},Object(m.b)("span",{parentName:"span",className:"pstrut",style:{height:"2.7em"}}),Object(m.b)("span",{parentName:"span",className:"sizing reset-size6 size3 mtight"},Object(m.b)("span",{parentName:"span",className:"mord mtight"},"2")))))))),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.16666666666666666em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"n"),Object(m.b)("span",{parentName:"span",className:"mclose"},")")))))," depth w.h.p. (using the\nimproved depth analysis of Fischer and Noever ","[2]","). More details about\nour algorithm implementation can be found in ","[3]","."),Object(m.b)("p",null,"The code for our implemenation is available\n",Object(m.b)("a",{parentName:"p",href:"https://github.com/ldhulipala/gbbs/tree/master/benchmarks/GraphColoring/Hasenplaugh14/"},"here"),"."),Object(m.b)("h2",{id:"cost-bounds"},"Cost Bounds"),Object(m.b)("p",null,"The algorithm runs in ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"O"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(m.b)("mi",{parentName:"mrow"},"m"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},")")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"O(m)")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.02778em"}},"O"),Object(m.b)("span",{parentName:"span",className:"mopen"},"("),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"m"),Object(m.b)("span",{parentName:"span",className:"mclose"},")")))))," work and ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"O"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(m.b)("msup",{parentName:"mrow"},Object(m.b)("mo",{parentName:"msup"},Object(m.b)("mi",{parentName:"mo"},"log"),Object(m.b)("mo",{parentName:"mo"},"\u2061")),Object(m.b)("mn",{parentName:"msup"},"2")),Object(m.b)("mi",{parentName:"mrow"},"n"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},")")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"O(\\log^2 n)")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1.148448em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.02778em"}},"O"),Object(m.b)("span",{parentName:"span",className:"mopen"},"("),Object(m.b)("span",{parentName:"span",className:"mop"},Object(m.b)("span",{parentName:"span",className:"mop"},"lo",Object(m.b)("span",{parentName:"span",style:{marginRight:"0.01389em"}},"g")),Object(m.b)("span",{parentName:"span",className:"msupsub"},Object(m.b)("span",{parentName:"span",className:"vlist-t"},Object(m.b)("span",{parentName:"span",className:"vlist-r"},Object(m.b)("span",{parentName:"span",className:"vlist",style:{height:"0.8984479999999999em"}},Object(m.b)("span",{parentName:"span",style:{top:"-3.1473400000000002em",marginRight:"0.05em"}},Object(m.b)("span",{parentName:"span",className:"pstrut",style:{height:"2.7em"}}),Object(m.b)("span",{parentName:"span",className:"sizing reset-size6 size3 mtight"},Object(m.b)("span",{parentName:"span",className:"mord mtight"},"2")))))))),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.16666666666666666em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"n"),Object(m.b)("span",{parentName:"span",className:"mclose"},")")))))," depth w.h.p."),Object(m.b)("h2",{id:"compiling-and-running"},"Compiling and Running"),Object(m.b)("p",null,"The benchmark can be compiled by running:"),Object(m.b)("pre",null,Object(m.b)("code",{parentName:"pre"},"bazel build -c opt //benchmarks/MaximalIndependentSet/RandomGreedy/...\n")),Object(m.b)("p",null,"It can then be run on a test input graph in the ",Object(m.b)("em",{parentName:"p"},"uncompressed format")," as follows:"),Object(m.b)("pre",null,Object(m.b)("code",{parentName:"pre"},"numactl -i all ./bazel-bin/benchmarks/MaximalIndependentSet/RandomGreedy/MaximalIndependentSet_main -s -m -src 1 inputs/rMatGraph_J_5_100\n")),Object(m.b)("p",null,"It can then be run on a test input graph in the ",Object(m.b)("em",{parentName:"p"},"compressed format")," as follows:"),Object(m.b)("pre",null,Object(m.b)("code",{parentName:"pre"},"numactl -i all ./bazel-bin/benchmarks/MaximalIndependentSet/RandomGreedy/MaximalIndependentSet_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda\n")),Object(m.b)("h2",{id:"references"},"References"),Object(m.b)("p",null,"[1]"," Guy Blelloch, Jeremy Fineman, and Julian Shun",Object(m.b)("br",null),"\n",Object(m.b)("em",{parentName:"p"},"Greedy Sequential Maximal Independent Set and Matching are Parallel on Average"),Object(m.b)("br",null),"\nProceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 308-317, 2012."),Object(m.b)("p",null,"[2]"," Manuella Fischer and Andreas Noever",Object(m.b)("br",null),"\n",Object(m.b)("em",{parentName:"p"},"Tight Analysis of Parallel Randomized Greedy MIS"),Object(m.b)("br",null),"\nProceedings of the ACM-SIAM Symposium on Discrete Algorithms (SODA), pp.2152-2160, 2018"),Object(m.b)("p",null,"[3]"," Laxman Dhulipala, Guy Blelloch, and Julian Shun",Object(m.b)("br",null),"\n",Object(m.b)("a",{parentName:"p",href:"https://ldhulipala.github.io/papers/gbbs_topc.pdf"},Object(m.b)("em",{parentName:"a"},"Theoretically Efficient Parallel Graph Algorithms Can Be Fast and Scalable")),Object(m.b)("br",null),"\nProceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 393-404, 2018."))}l.isMDXComponent=!0}}]);