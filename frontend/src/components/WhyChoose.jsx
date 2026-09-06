import React, { useEffect, useRef, useState } from 'react';
import { ArrowRight, Sparkles, Database, Layers } from 'lucide-react';
import './whyChoose.css';

const features = [
  {
    icon: Sparkles,
    title: "Turn raw eDNA into actionable biodiversity insights",
    points: [
      "Collect, process, and interpret environmental DNA samples in a single streamlined platform.",
      "Eliminate guesswork with automated quality control, denoising, and contamination detection.",
      "Leverage AI-driven taxonomic classification to reveal hidden species."
    ],
    linkText: "Read more",
    videoSrc: "/assets/videos/edna-insights-demo.mp4",
    videoPoster: "/assets/videos/edna-insights-poster.jpg"
  },
  {
    icon: Database,
    title: "Build a trusted biodiversity knowledge base",
    points: [
      "Track taxonomic assignments, community structures, and diversity indices across projects.",
      "Design unified dashboards that connect sequencing results, metadata, and annotations.",
      "Empower researchers with reproducible, shareable outputs for downstream science."
    ],
    linkText: "Read more",
    videoSrc: "/assets/videos/knowledge-base-demo.mp4",
    videoPoster: "/assets/videos/knowledge-base-poster.jpg"
  },
  {
    icon: Layers,
    title: "Evolve with your environmental research needs",
    points: [
      "Access modular workflows for QC, novelty detection, clustering, and phylogeny.",
      "Integrate seamlessly with existing sequencing pipelines, databases, and visualization tools.",
      "Scale on a secure, flexible platform designed to grow with your datasets and discoveries."
    ],
    linkText: "Read more",
    videoSrc: "/assets/videos/evolve-platform-demo.mp4",
    videoPoster: "/assets/videos/evolve-platform-poster.jpg"
  }
];

export default function WhyChoose() {
  const [visibleCards, setVisibleCards] = useState([]);
  const [activeVideo, setActiveVideo] = useState(0);
  const cardRefs = useRef([]);
  const sectionRef = useRef(null);
  const videoRefs = useRef([]);
  const cardVideoRefs = useRef([]); // New ref for individual card videos

  useEffect(() => {
    const observer = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => {
          if (entry.isIntersecting) {
            const index = cardRefs.current.indexOf(entry.target);
            if (index !== -1 && !visibleCards.includes(index)) {
              setTimeout(() => {
                setVisibleCards(prev => [...prev, index]);
              }, index * 150); // Stagger animation
            }
          }
        });
      },
      { threshold: 0.1, rootMargin: '0px 0px -50px 0px' }
    );

    cardRefs.current.forEach((ref) => {
      if (ref) observer.observe(ref);
    });

    return () => observer.disconnect();
  }, [visibleCards]);

  const handleCardHover = (index) => {
    setActiveVideo(index);
    // Play the video in the card on desktop
    if (window.innerWidth >= 969 && cardVideoRefs.current[index]) {
      cardVideoRefs.current[index].play();
    }
    // Play the video in the stack on mobile
    if (window.innerWidth < 969 && videoRefs.current[index]) {
      videoRefs.current[index].play();
    }
  };

  const handleCardLeave = (index) => {
    // Pause the video when hover leaves
    if (cardVideoRefs.current[index]) {
      cardVideoRefs.current[index].pause();
    }
    if (videoRefs.current[index]) {
      videoRefs.current[index].pause();
    }
  };

  return (
    <section className="why-choose-section" ref={sectionRef}>
      <div className="why-choose-container">
        <div className="why-choose-header">
          <h2 className="why-choose-title">Why choose AQUADEX</h2>
          <p className="why-choose-subtitle">
            Transform your environmental DNA research with our comprehensive platform
          </p>
        </div>

        <div className="why-choose-content">
          {/* Features Column */}
          <div className="features-column">
            {features.map((feature, index) => {
              const Icon = feature.icon;
              return (
                <div
                  key={index}
                  ref={el => cardRefs.current[index] = el}
                  className={`feature-card ${visibleCards.includes(index) ? 'fade-up' : ''} ${activeVideo === index ? 'active' : ''}`}
                  onMouseEnter={() => handleCardHover(index)}
                  onMouseLeave={() => handleCardLeave(index)}
                  onClick={() => setActiveVideo(index)}
                >
                  <div className="feature-content">
                    <div className="feature-icon">
                      <Icon size={24} />
                    </div>
                    <h3 className="feature-title">{feature.title}</h3>
                    <ul className="feature-points">
                      {feature.points.map((point, pointIndex) => (
                        <li key={pointIndex} className="feature-point">
                          {point}
                        </li>
                      ))}
                    </ul>
                    <a href="#" className="feature-link">
                      {feature.linkText}
                      <ArrowRight size={16} className="link-arrow" />
                    </a>
                  </div>
                  
                  {/* Individual video container for desktop */}
                  <div className="feature-video-container">
                    <video
                      ref={el => cardVideoRefs.current[index] = el}
                      className="feature-video"
                      muted
                      loop
                      playsInline
                      poster={feature.videoPoster}
                      aria-label={`Demo: ${feature.title}`}
                    >
                      <source src={feature.videoSrc} type="video/mp4" />
                      Your browser does not support the video tag.
                    </video>
                  </div>
                </div>
              );
            })}
          </div>

          {/* Video Column - Only visible on mobile */}
          <div className="video-column">
            <div className="video-stack">
              {features.map((feature, index) => (
                <div
                  key={index}
                  className={`video-wrapper ${activeVideo === index ? 'active' : ''}`}
                >
                  <video
                    ref={el => videoRefs.current[index] = el}
                    className="demo-video"
                    muted
                    loop
                    playsInline
                    poster={feature.videoPoster}
                    aria-label={`Demo: ${feature.title}`}
                  >
                    <source src={feature.videoSrc} type="video/mp4" />
                    Your browser does not support the video tag.
                  </video>
                  <div className="video-caption">
                    {feature.title}
                  </div>
                </div>
              ))}
            </div>
          </div>
        </div>
      </div>
    </section>
  );
}